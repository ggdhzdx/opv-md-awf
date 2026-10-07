"""
bulkgraph_core.geometry: SVD plane fitting, vector geometry, and periodic boundary helpers.
"""

from typing import Tuple, Optional
import numpy as np


def calc_plane_normal(coords: np.ndarray) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Fits a best-fit plane to a set of heavy atom coordinates using SVD.

    Parameters
    ----------
    coords : np.ndarray, shape (N, 3)
        Cartesian coordinates in nanometers.

    Returns
    -------
    normal : np.ndarray, shape (3,)
        Unit normal vector of the fitted plane.
    centroid : np.ndarray, shape (3,)
        Geometric center of the input coordinates.
    rmsd : float
        Out-of-plane RMSD in nanometers.
    """
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError(f"coords must have shape (N, 3), got {coords.shape}")
    n_atoms = coords.shape[0]
    if n_atoms < 3:
        raise ValueError(f"At least 3 atoms required for plane fitting, got {n_atoms}")

    centroid = np.mean(coords, axis=0)
    centered = coords - centroid

    # SVD: centered = U @ diag(s) @ Vt
    # The normal vector corresponds to the right singular vector with smallest singular value
    _, s, vt = np.linalg.svd(centered, full_matrices=False)
    normal = vt[2, :].copy()

    norm_val = np.linalg.norm(normal)
    if norm_val < 1e-12:
        raise ValueError("Collinear or degenerate coordinates; normal vector norm is zero.")
    normal = normal / norm_val

    rmsd = float(s[2] / np.sqrt(n_atoms))
    return normal, centroid, rmsd


def calc_normal_angle(n1: np.ndarray, n2: np.ndarray) -> float:
    """
    Calculates the undirected tilt angle in degrees between two plane normal vectors.
    Since planes have no intrinsic directionality, cos(theta) = |n1 . n2|.

    Returns
    -------
    angle_deg : float in [0.0, 90.0] degrees.
    """
    dot = float(np.dot(n1, n2))
    cos_theta = np.clip(np.abs(dot), 0.0, 1.0)
    angle_rad = np.arccos(cos_theta)
    return float(np.degrees(angle_rad))


def align_normal_signs(n1: np.ndarray, n2: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Sign-aligns n2 with n1 such that n1 . n2_aligned >= 0, preventing
    anti-parallel vector cancellation when taking their mean.

    Returns
    -------
    n2_aligned : np.ndarray, shape (3,)
        Sign-aligned n2.
    n_avg : np.ndarray, shape (3,)
        Normalized average unit normal with canonical sign.
    """
    dot = float(np.dot(n1, n2))
    if dot < 0.0:
        n2_aligned = -n2
    else:
        n2_aligned = n2

    n_sum = n1 + n2_aligned
    norm_sum = np.linalg.norm(n_sum)
    if norm_sum < 1e-12:
        raise ValueError("Normal vector cancellation detected after sign alignment.")
    n_avg = n_sum / norm_sum

    # Canonicalize direction: first non-negligible component > 0
    for comp in n_avg:
        if np.abs(comp) > 1e-8:
            if comp < 0:
                n_avg = -n_avg
            break

    return n2_aligned, n_avg


def pbc_displacement(r1: np.ndarray, r2: np.ndarray, box: np.ndarray) -> np.ndarray:
    """
    Computes displacement vector (r1 - r2) under periodic boundary conditions
    for an orthorhombic simulation box.

    Parameters
    ----------
    r1, r2 : np.ndarray, shape (3,) or (N, 3)
    box : np.ndarray
        Box dimensions [Lx, Ly, Lz] or [Lx, Ly, Lz, alpha, beta, gamma].

    Returns
    -------
    disp_pbc : np.ndarray
    """
    box_lengths = np.asarray(box[:3], dtype=np.float64)
    dr = r1 - r2
    dr_pbc = dr - box_lengths * np.round(dr / box_lengths)
    return dr_pbc


def calc_normal_projection(
    c1: np.ndarray,
    c2: np.ndarray,
    n_avg: np.ndarray,
    box: Optional[np.ndarray] = None
) -> float:
    """
    Calculates the projection distance between two centroids along the mean plane normal.

    Parameters
    ----------
    c1, c2 : np.ndarray, shape (3,)
        Centroid coordinates in nanometers.
    n_avg : np.ndarray, shape (3,)
        Mean unit normal vector.
    box : Optional[np.ndarray]
        Box dimensions for PBC minimum image displacement.

    Returns
    -------
    proj : float
        Absolute projection distance in nanometers.
    """
    if box is not None:
        dr = pbc_displacement(c1, c2, box)
    else:
        dr = c1 - c2
    proj = float(np.abs(np.dot(n_avg, dr)))
    return proj


def pbc_distance_matrix(
    coords1: np.ndarray,
    coords2: np.ndarray,
    box: np.ndarray
) -> np.ndarray:
    """
    Computes pairwise Euclidean distance matrix under PBC minimum image convention.

    Parameters
    ----------
    coords1 : np.ndarray, shape (N, 3)
    coords2 : np.ndarray, shape (M, 3)
    box : np.ndarray
        Box dimensions [Lx, Ly, Lz, 90, 90, 90] or [Lx, Ly, Lz].

    Returns
    -------
    dists : np.ndarray, shape (N, M)
    """
    box_array = np.asarray(box, dtype=np.float64)
    if box_array.shape not in {(3,), (6,)}:
        raise ValueError(
            "BulkGraph currently requires an orthorhombic box as "
            "[Lx,Ly,Lz] or [Lx,Ly,Lz,90,90,90]"
        )
    if box_array.shape == (6,) and not np.allclose(box_array[3:], 90.0, atol=1e-6):
        raise ValueError("BulkGraph currently requires orthorhombic box angles")
    lengths = box_array[:3]
    if np.any(lengths <= 0):
        raise ValueError("box lengths must be positive")
    c1 = np.asarray(coords1, dtype=np.float64)
    c2 = np.asarray(coords2, dtype=np.float64)
    delta = c1[:, None, :] - c2[None, :, :]
    delta -= lengths * np.round(delta / lengths)
    return np.linalg.norm(delta, axis=-1)
