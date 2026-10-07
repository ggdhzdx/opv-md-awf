"""
bulkgraph_core.stacking: Multi-descriptor geometric stacking evaluation and pure NumPy grid overlap area.
"""

from typing import Dict, Any, List, Optional
import numpy as np
from .geometry import (
    calc_plane_normal,
    calc_normal_angle,
    align_normal_signs,
    calc_normal_projection,
    pbc_displacement,
    pbc_distance_matrix,
)

DEFAULT_RADII_NM = {
    "C": 0.170,
    "N": 0.155,
    "O": 0.152,
    "F": 0.147,
    "S": 0.180,
}


class StackingGateError(Exception):
    """Raised when stacking area evaluation is requested without satisfying the stacking gate."""
    pass


class BlockedComponentError(Exception):
    """Raised when an unbuilt or blocked engine is called."""
    pass


def evaluate_stacking_pair(
    coords_A: np.ndarray,
    coords_B: np.ndarray,
    box: np.ndarray,
    dist_thresh: float = 0.36,
    angle_thresh: float = 30.0,
    proj_thresh: float = 0.25,
) -> Dict[str, Any]:
    """
    Evaluates geometric stacking criteria for a pair of conjugated fragments.

    Parameters
    ----------
    coords_A, coords_B : np.ndarray, shape (N, 3) and (M, 3)
        Heavy-atom coordinates in nanometers.
    box : np.ndarray
        Simulation box dimensions.
    dist_thresh : float
        Primary heavy-atom distance threshold in nanometers (default 0.36 nm).
    angle_thresh : float
        Planar normal tilt angle threshold in degrees (default 30.0 deg).
    proj_thresh : float
        Normal projection distance threshold in nanometers (default 0.25 nm).

    Returns
    -------
    result : dict
        Contains 'is_stacked', 'is_stacked_sens_040', 'angle_deg', 'proj_nm',
        'min_dist_nm', 'normal_A', 'normal_B', 'normal_avg', 'centroid_A', 'centroid_B'.
    """
    normal_A, centroid_A, _ = calc_plane_normal(coords_A)
    normal_B, centroid_B, _ = calc_plane_normal(coords_B)

    angle_deg = calc_normal_angle(normal_A, normal_B)
    _, normal_avg = align_normal_signs(normal_A, normal_B)

    proj_nm = calc_normal_projection(centroid_A, centroid_B, normal_avg, box=box)

    dists = pbc_distance_matrix(coords_A, coords_B, box)
    min_dist_nm = float(np.min(dists))

    angle_ok = (angle_deg <= angle_thresh)
    proj_ok = (proj_nm >= proj_thresh)
    dist_primary_ok = (min_dist_nm <= dist_thresh)
    dist_sens_ok = (min_dist_nm <= 0.40)

    is_stacked = bool(angle_ok and proj_ok and dist_primary_ok)
    is_stacked_sens_040 = bool(angle_ok and proj_ok and dist_sens_ok)

    return {
        "is_stacked": is_stacked,
        "is_stacked_sens_040": is_stacked_sens_040,
        "angle_deg": float(angle_deg),
        "proj_nm": float(proj_nm),
        "min_dist_nm": float(min_dist_nm),
        "normal_A": normal_A,
        "normal_B": normal_B,
        "normal_avg": normal_avg,
        "centroid_A": centroid_A,
        "centroid_B": centroid_B,
    }


def projected_cofacial_overlap_area(
    coords_A: np.ndarray,
    elements_A: List[str],
    coords_B: np.ndarray,
    elements_B: List[str],
    box: Optional[np.ndarray] = None,
    gate_passed: bool = False,
    h_init: float = 0.01,
    tol_abs: float = 0.001,
    tol_rel: float = 0.002,
    radii_table: Optional[Dict[str, float]] = None,
) -> Dict[str, Any]:
    """
    Computes the projected cofacial overlap area between two stacked fragments
    using deterministic pure NumPy cell-centered grid integration with adaptive refinement.

    Parameters
    ----------
    coords_A, coords_B : np.ndarray, shape (N, 3), (M, 3)
        Heavy-atom coordinates in nanometers.
    elements_A, elements_B : list of str
        Element symbols for heavy atoms.
    box : Optional[np.ndarray]
        Simulation box dimensions for PBC coordinate alignment.
    gate_passed : bool
        Must be True (from geometric stacking filter). If False, raises StackingGateError.
    h_init : float
        Initial grid spacing in nanometers (default 0.01 nm).
    tol_abs : float
        Absolute tolerance in nm^2 (default 0.001 nm^2).
    tol_rel : float
        Relative tolerance (default 0.002).
    radii_table : Optional[Dict[str, float]]
        van der Waals radii mapping in nanometers.

    Returns
    -------
    dict with:
        'status': 'CONVERGED' or 'BLOCKED_NUMERICAL_CONVERGENCE'
        'area_nm2': float or None
        'refinement_level': int
        'h_final_nm': float
    """
    if not gate_passed:
        raise StackingGateError(
            "Stacking gate contract failed: gate_passed must be True to evaluate projected cofacial overlap area."
        )

    # Element validation & hydrogen policy
    if "H" in elements_A or "H" in elements_B:
        raise ValueError("Hydrogen atoms must be strictly excluded from projected cofacial overlap area.")

    radii = radii_table or DEFAULT_RADII_NM
    try:
        r_A = np.array([radii[el] for el in elements_A], dtype=np.float64)
        r_B = np.array([radii[el] for el in elements_B], dtype=np.float64)
    except KeyError as e:
        raise KeyError(f"Unrecognized element in radii table: {e}. Fail closed.")

    # Plane normals and centroids
    normal_A, centroid_A, _ = calc_plane_normal(coords_A)
    normal_B, centroid_B, _ = calc_plane_normal(coords_B)
    _, normal_avg = align_normal_signs(normal_A, normal_B)

    # Place fragment B in the unwrapped local frame of fragment A across PBC
    if box is not None:
        dr_c = pbc_displacement(centroid_B, centroid_A, box)
        coords_B_local = coords_B + (centroid_A + dr_c - centroid_B)
        centroid_B_local = centroid_A + dr_c
    else:
        coords_B_local = coords_B.copy()
        centroid_B_local = centroid_B

    c0 = 0.5 * (centroid_A + centroid_B_local)

    # Construct orthonormal 2D basis (u1, u2) on mean plane
    # Select cardinal axis of minimum normal component
    cardinal = np.zeros(3)
    cardinal[np.argmin(np.abs(normal_avg))] = 1.0
    v1 = cardinal - np.dot(cardinal, normal_avg) * normal_avg
    u1 = v1 / np.linalg.norm(v1)
    u2 = np.cross(normal_avg, u1)
    u2 = u2 / np.linalg.norm(u2)

    # Project coordinates to 2D
    rel_A = coords_A - c0
    xi_A = np.dot(rel_A, u1)
    eta_A = np.dot(rel_A, u2)

    rel_B = coords_B_local - c0
    xi_B = np.dot(rel_B, u1)
    eta_B = np.dot(rel_B, u2)

    # 2D Bounding boxes
    xmin_A, xmax_A = np.min(xi_A - r_A), np.max(xi_A + r_A)
    ymin_A, ymax_A = np.min(eta_A - r_A), np.max(eta_A + r_A)

    xmin_B, xmax_B = np.min(xi_B - r_B), np.max(xi_B + r_B)
    ymin_B, ymax_B = np.min(eta_B - r_B), np.max(eta_B + r_B)

    # Overlap bounding box
    x_min = max(xmin_A, xmin_B)
    x_max = min(xmax_A, xmax_B)
    y_min = max(ymin_A, ymin_B)
    y_max = min(ymax_A, ymax_B)

    if x_max <= x_min or y_max <= y_min:
        return {
            "status": "CONVERGED",
            "area_nm2": 0.0,
            "refinement_level": 0,
            "h_final_nm": h_init,
        }

    # Helper function for grid integration at spacing h
    r_A_sq = r_A ** 2
    r_B_sq = r_B ** 2

    def _grid_area(h: float) -> float:
        xs = np.arange(x_min + 0.5 * h, x_max, h)
        ys = np.arange(y_min + 0.5 * h, y_max, h)
        if len(xs) == 0 or len(ys) == 0:
            return 0.0

        # Broadcast grid points (ny, nx)
        grid_x = xs[None, :]  # shape (1, nx)
        grid_y = ys[:, None]  # shape (ny, 1)

        # Check membership in union A: distance_sq <= r_A^2 for any atom i
        # To be memory efficient, evaluate chunked if needed, but overlap box is small (~0.5 - 2 nm)
        # Typically nx, ny ~ 50-200 points
        in_A = np.zeros((len(ys), len(xs)), dtype=bool)
        for i in range(len(xi_A)):
            dx = grid_x - xi_A[i]
            dy = grid_y - eta_A[i]
            in_A |= ((dx * dx + dy * dy) <= r_A_sq[i])

        in_B = np.zeros((len(ys), len(xs)), dtype=bool)
        for j in range(len(xi_B)):
            dx = grid_x - xi_B[j]
            dy = grid_y - eta_B[j]
            in_B |= ((dx * dx + dy * dy) <= r_B_sq[j])

        overlap_cells = np.count_nonzero(in_A & in_B)
        return float(overlap_cells * (h * h))

    # Refinement schedule: h -> h/2 -> h/4
    a_h = _grid_area(h_init)
    h_mid = h_init / 2.0
    a_h2 = _grid_area(h_mid)

    tol1 = max(tol_abs, tol_rel * a_h2)
    if abs(a_h2 - a_h) <= tol1:
        return {
            "status": "CONVERGED",
            "area_nm2": float(a_h2),
            "refinement_level": 1,
            "h_final_nm": float(h_mid),
        }

    h_fine = h_init / 4.0
    a_h4 = _grid_area(h_fine)
    tol2 = max(tol_abs, tol_rel * a_h4)
    if abs(a_h4 - a_h2) <= tol2:
        return {
            "status": "CONVERGED",
            "area_nm2": float(a_h4),
            "refinement_level": 2,
            "h_final_nm": float(h_fine),
        }

    # Unconverged after h/4: fail closed
    return {
        "status": "BLOCKED_NUMERICAL_CONVERGENCE",
        "area_nm2": None,
        "refinement_level": 2,
        "h_final_nm": float(h_fine),
    }


def exact_analytical_overlap_area(*args, **kwargs) -> Any:
    """
    Exact analytical area engine is unbuilt and remains strictly blocked.
    """
    raise BlockedComponentError(
        "Exact analytical area engine is unbuilt and remains BLOCKED_PENDING_EXACT_AREA_ENGINE."
    )
