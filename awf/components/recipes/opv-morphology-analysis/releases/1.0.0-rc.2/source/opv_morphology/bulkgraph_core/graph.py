"""
bulkgraph_core.graph: Dual-graph network generators (G_contact and G_stack) and stacking motif evaluation.
"""

from typing import Dict, Any, List, Tuple, Optional
import numpy as np
from .geometry import (
    calc_plane_normal,
    calc_normal_angle,
    align_normal_signs,
    calc_normal_projection,
    pbc_displacement,
    pbc_distance_matrix,
)
from .stacking import evaluate_stacking_pair, projected_cofacial_overlap_area
from .topology import compute_graph_metrics


def build_contact_graph(
    mol_conjugated_coords: List[np.ndarray],
    box: np.ndarray,
    cutoff_nm: float = 0.45,
    sens_cutoffs_nm: Optional[List[float]] = None,
) -> Dict[str, Any]:
    """
    Constructs the molecular contact graph G_contact among NFA molecules based on
    conjugated heavy-atom minimum PBC distance.

    Parameters
    ----------
    mol_conjugated_coords : list of np.ndarray
        List of length N_mols (e.g. 264), each array of shape (N_conj_heavy, 3).
    box : np.ndarray
        Simulation box [Lx, Ly, Lz, 90, 90, 90] or [Lx, Ly, Lz].
    cutoff_nm : float
        Primary contact distance cutoff in nanometers (default 0.45 nm).
    sens_cutoffs_nm : Optional[List[float]]
        Sensitivity cutoffs (e.g. [0.40, 0.50]).

    Returns
    -------
    dict with primary and sensitivity graph metrics and edges.
    """
    n_mols = len(mol_conjugated_coords)
    sens_cutoffs = sens_cutoffs_nm or [0.40, 0.50]

    # Pre-calculate molecule conjugated centroids and approximate bounding radii for filtering
    mol_centroids = np.array([np.mean(coords, axis=0) for coords in mol_conjugated_coords], dtype=np.float64)
    mol_radii = np.array([np.max(np.linalg.norm(coords - c, axis=1)) for coords, c in zip(mol_conjugated_coords, mol_centroids)], dtype=np.float64)

    edges_primary: List[Tuple[int, int]] = []
    edges_sens: Dict[float, List[Tuple[int, int]]] = {sc: [] for sc in sens_cutoffs}
    max_sens = max(max(sens_cutoffs), cutoff_nm)

    box_l = box[:3]
    diff = mol_centroids[:, np.newaxis, :] - mol_centroids[np.newaxis, :, :]
    diff -= np.round(diff / box_l) * box_l
    dist_mat = np.linalg.norm(diff, axis=-1)

    threshold_mat = mol_radii[:, np.newaxis] + mol_radii[np.newaxis, :] + max_sens
    u_cand, v_cand = np.where(np.triu(dist_mat <= threshold_mat, k=1))

    for u, v in zip(u_cand, v_cand):
        coords_u = mol_conjugated_coords[u]
        dists = pbc_distance_matrix(coords_u, mol_conjugated_coords[v], box)
        min_d = float(np.min(dists))

        if min_d <= cutoff_nm:
            edges_primary.append((int(u), int(v)))
        for sc in sens_cutoffs:
            if min_d <= sc:
                edges_sens[sc].append((int(u), int(v)))

    primary_metrics = compute_graph_metrics(n_mols, edges_primary)
    result = {
        "cutoff_nm": cutoff_nm,
        "edge_count": len(edges_primary),
        "edges": edges_primary,
        "metrics": primary_metrics,
        "sensitivities": {},
    }

    for sc in sens_cutoffs:
        result["sensitivities"][str(sc)] = {
            "cutoff_nm": sc,
            "edge_count": len(edges_sens[sc]),
            "metrics": compute_graph_metrics(n_mols, edges_sens[sc]),
        }

    return result


def build_stacking_graph(
    mol_fragment_dict_list: List[Dict[str, Dict[str, Any]]],
    box: np.ndarray,
    dist_thresh: float = 0.36,
    angle_thresh: float = 30.0,
    proj_thresh: float = 0.25,
    compute_overlap_area: bool = True,
) -> Dict[str, Any]:
    """
    Constructs the BulkGraph stacking graph G_stack and evaluates motif-specific
    stacking counts (IC-IC, IC-core, core-core) and participation fraction.

    Parameters
    ----------
    mol_fragment_dict_list : list of dict
        Length N_mols (e.g. 264). Each entry is a dict:
        {
            'ADN': {'coords': np.ndarray (27, 3), 'elements': list of str (len 27)},
            'AAC_left': {'coords': np.ndarray (18, 3), 'elements': list of str (len 18)},
            'AAC_right': {'coords': np.ndarray (18, 3), 'elements': list of str (len 18)},
        }
    box : np.ndarray
        Simulation box.
    dist_thresh : float
        Primary stacking distance (default 0.36 nm).
    angle_thresh : float
        Tilt angle threshold (default 30.0 deg).
    proj_thresh : float
        Normal projection distance threshold (default 0.25 nm).
    compute_overlap_area : bool
        Whether to compute projected cofacial overlap area on primary stacked pairs.

    Returns
    -------
    dict with stacking edges, motif counts, f_atom_stack, and topological metrics.
    """
    n_mols = len(mol_fragment_dict_list)

    # Precompute planar SVD for each fragment on each molecule
    frag_names = ["ADN", "AAC_left", "AAC_right"]
    mol_frag_geometry = []
    total_conj_heavy_atoms = 0

    for u, m_frags in enumerate(mol_fragment_dict_list):
        u_geom = {}
        for fname in frag_names:
            coords = m_frags[fname]["coords"]
            elements = m_frags[fname]["elements"]
            normal, centroid, rmsd = calc_plane_normal(coords)
            u_geom[fname] = {
                "coords": coords,
                "elements": elements,
                "normal": normal,
                "centroid": centroid,
                "rmsd": rmsd,
            }
            if u == 0:
                pass  # counting total per mol later
        mol_frag_geometry.append(u_geom)

    # Count total conjugated heavy atoms across all molecules
    total_conj_heavy_atoms = sum(
        sum(len(mol_frag_geometry[u][fname]["coords"]) for fname in frag_names)
        for u in range(n_mols)
    )

    edges_primary: List[Tuple[int, int]] = []
    edges_sens_040: List[Tuple[int, int]] = []
    edges_angle_025: List[Tuple[int, int]] = []
    edges_angle_035: List[Tuple[int, int]] = []

    # Motif counts
    motif_counts_primary = {
        "IC_IC": 0,
        "IC_core": 0,
        "core_core": 0,
    }
    motif_counts_sens_040 = {
        "IC_IC": 0,
        "IC_core": 0,
        "core_core": 0,
    }
    motif_counts_angle_025 = {
        "IC_IC": 0,
        "IC_core": 0,
        "core_core": 0,
    }
    motif_counts_angle_035 = {
        "IC_IC": 0,
        "IC_core": 0,
        "core_core": 0,
    }

    # Stacking overlap areas per motif
    motif_areas: Dict[str, List[float]] = {
        "IC_IC": [],
        "IC_core": [],
        "core_core": [],
    }

    # Precompute molecular conjugated centroids and bounding radii for triangle-inequality filtering
    mol_centroids = []
    mol_radii = []
    for u in range(n_mols):
        all_u_coords = np.vstack([mol_frag_geometry[u][fname]["coords"] for fname in frag_names])
        c_u = np.mean(all_u_coords, axis=0)
        r_u = float(np.max(np.linalg.norm(all_u_coords - c_u, axis=1)))
        mol_centroids.append(c_u)
        mol_radii.append(r_u)
    mol_centroids = np.array(mol_centroids, dtype=np.float64)
    mol_radii = np.array(mol_radii, dtype=np.float64)

    # Stacking participation tracking: set of (mol_idx, fname, atom_local_idx)
    stacking_participant_atoms = set()

    box_l = box[:3]
    diff = mol_centroids[:, np.newaxis, :] - mol_centroids[np.newaxis, :]
    diff -= np.round(diff / box_l) * box_l
    dist_mat = np.linalg.norm(diff, axis=-1)

    threshold_mat = mol_radii[:, np.newaxis] + mol_radii[np.newaxis, :] + 0.40
    u_cand, v_cand = np.where(np.triu(dist_mat <= threshold_mat, k=1))

    for u, v in zip(u_cand, v_cand):
        u = int(u)
        v = int(v)
        u_geom = mol_frag_geometry[u]
        v_geom = mol_frag_geometry[v]
        mol_stacked_primary = False
        mol_stacked_sens = False
        mol_stacked_angle_025 = False
        mol_stacked_angle_035 = False

        for f_u_name in frag_names:
            fu = u_geom[f_u_name]
            for f_v_name in frag_names:
                fv = v_geom[f_v_name]

                # 1. Tilt angle candidate filter (evaluate when theta <= 35 deg)
                angle = calc_normal_angle(fu["normal"], fv["normal"])
                if angle > 35.0:
                    continue

                # 2. Normal projection filter with sign alignment
                _, n_avg = align_normal_signs(fu["normal"], fv["normal"])
                proj = calc_normal_projection(fu["centroid"], fv["centroid"], n_avg, box=box)
                if proj < proj_thresh:
                    continue

                # 3. Minimum heavy-atom distance under PBC
                dists = pbc_distance_matrix(fu["coords"], fv["coords"], box)
                min_dist = float(np.min(dists))

                # Classify motif type
                if f_u_name == "ADN" and f_v_name == "ADN":
                    m_type = "core_core"
                elif f_u_name != "ADN" and f_v_name != "ADN":
                    m_type = "IC_IC"
                else:
                    m_type = "IC_core"

                # angle_025: theta <= 25 and proj >= 0.25 and d_min <= 0.36
                if angle <= 25.0 and min_dist <= dist_thresh:
                    mol_stacked_angle_025 = True
                    motif_counts_angle_025[m_type] += 1

                # Primary stacking evaluation: theta <= 30 and proj >= 0.25 and d_min <= 0.36
                if angle <= angle_thresh and min_dist <= dist_thresh:
                    mol_stacked_primary = True
                    motif_counts_primary[m_type] += 1

                    # Track participation
                    close_u = np.where(np.min(dists, axis=1) <= dist_thresh)[0]
                    for idx_u in close_u:
                        stacking_participant_atoms.add((u, f_u_name, idx_u))
                    close_v = np.where(np.min(dists, axis=0) <= dist_thresh)[0]
                    for idx_v in close_v:
                        stacking_participant_atoms.add((v, f_v_name, idx_v))

                    # Projected cofacial overlap area
                    if compute_overlap_area:
                        area_res = projected_cofacial_overlap_area(
                            fu["coords"],
                            fu["elements"],
                            fv["coords"],
                            fv["elements"],
                            box=box,
                            gate_passed=True,
                        )
                        if area_res["status"] == "CONVERGED" and area_res["area_nm2"] is not None:
                            motif_areas[m_type].append(area_res["area_nm2"])

                # angle_035: theta <= 35 and proj >= 0.25 and d_min <= 0.36
                if angle <= 35.0 and min_dist <= dist_thresh:
                    mol_stacked_angle_035 = True
                    motif_counts_angle_035[m_type] += 1

                # Sensitivity stacking evaluation: theta <= 30 and proj >= 0.25 and d_min <= 0.40
                if angle <= angle_thresh and min_dist <= 0.40:
                    mol_stacked_sens = True
                    motif_counts_sens_040[m_type] += 1

        if mol_stacked_angle_025:
            edges_angle_025.append((u, v))
        if mol_stacked_primary:
            edges_primary.append((u, v))
        if mol_stacked_angle_035:
            edges_angle_035.append((u, v))
        if mol_stacked_sens:
            edges_sens_040.append((u, v))

    primary_metrics = compute_graph_metrics(n_mols, edges_primary)
    sens_metrics = compute_graph_metrics(n_mols, edges_sens_040)
    angle_025_metrics = compute_graph_metrics(n_mols, edges_angle_025)
    angle_035_metrics = compute_graph_metrics(n_mols, edges_angle_035)

    f_atom_stack = float(len(stacking_participant_atoms) / total_conj_heavy_atoms) if total_conj_heavy_atoms > 0 else 0.0

    mean_areas = {
        m_type: float(np.mean(motif_areas[m_type])) if len(motif_areas[m_type]) > 0 else 0.0
        for m_type in motif_areas
    }

    return {
        "primary": {
            "dist_thresh_nm": dist_thresh,
            "angle_thresh_deg": angle_thresh,
            "proj_thresh_nm": proj_thresh,
            "edge_count": len(edges_primary),
            "edges": edges_primary,
            "metrics": primary_metrics,
            "motifs": motif_counts_primary,
            "f_atom_stack": f_atom_stack,
            "mean_overlap_area_nm2": mean_areas,
            "overlap_area_counts": {m: len(motif_areas[m]) for m in motif_areas},
        },
        "sensitivity_040": {
            "dist_thresh_nm": 0.40,
            "edge_count": len(edges_sens_040),
            "edges": edges_sens_040,
            "metrics": sens_metrics,
            "motifs": motif_counts_sens_040,
        },
        "angle_025": {
            "dist_thresh_nm": dist_thresh,
            "angle_thresh_deg": 25.0,
            "proj_thresh_nm": proj_thresh,
            "edge_count": len(edges_angle_025),
            "edges": edges_angle_025,
            "metrics": angle_025_metrics,
            "motifs": motif_counts_angle_025,
        },
        "angle_035": {
            "dist_thresh_nm": dist_thresh,
            "angle_thresh_deg": 35.0,
            "proj_thresh_nm": proj_thresh,
            "edge_count": len(edges_angle_035),
            "edges": edges_angle_035,
            "metrics": angle_035_metrics,
            "motifs": motif_counts_angle_035,
        },
    }
