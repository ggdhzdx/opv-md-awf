"""
bulkgraph_core: Trajectory-agnostic, chemistry-agnostic geometric and topological library.
"""

from .geometry import (
    calc_plane_normal,
    calc_normal_angle,
    align_normal_signs,
    pbc_displacement,
    calc_normal_projection,
    pbc_distance_matrix,
)
from .stacking import (
    evaluate_stacking_pair,
    projected_cofacial_overlap_area,
    StackingGateError,
    BlockedComponentError,
)
from .graph import build_contact_graph, build_stacking_graph
from .topology import compute_graph_metrics

__all__ = [
    "calc_plane_normal",
    "calc_normal_angle",
    "align_normal_signs",
    "pbc_displacement",
    "calc_normal_projection",
    "pbc_distance_matrix",
    "evaluate_stacking_pair",
    "projected_cofacial_overlap_area",
    "StackingGateError",
    "BlockedComponentError",
    "build_contact_graph",
    "build_stacking_graph",
    "compute_graph_metrics",
]
