"""
bulkgraph_core.topology: Pure topological graph algorithms (degree sequences, BFS connected components, LCC, isolated fractions).
"""

from typing import Dict, Any, List, Set, Union, Tuple
from collections import deque
import numpy as np


def compute_connected_components_bfs(n_nodes: int, adj_list: List[List[int]]) -> List[List[int]]:
    """
    Computes connected components of an undirected graph using Breadth-First Search (BFS).

    Parameters
    ----------
    n_nodes : int
        Total number of nodes (0..n_nodes-1).
    adj_list : List[List[int]]
        Adjacency list where adj_list[u] contains neighbors of node u.

    Returns
    -------
    components : List[List[int]]
        List of connected components (each a list of node indices), sorted by size descending.
    """
    visited = [False] * n_nodes
    components = []

    for start_node in range(n_nodes):
        if not visited[start_node]:
            comp = []
            queue = deque([start_node])
            visited[start_node] = True
            while queue:
                u = queue.popleft()
                comp.append(u)
                for v in adj_list[u]:
                    if not visited[v]:
                        visited[v] = True
                        queue.append(v)
            components.append(comp)

    components.sort(key=lambda c: len(c), reverse=True)
    return components


def compute_graph_metrics(
    n_nodes: int,
    edges_or_adj: Union[List[Tuple[int, int]], np.ndarray, List[List[int]]]
) -> Dict[str, Any]:
    """
    Computes standard topological descriptors for an undirected molecular graph.

    Parameters
    ----------
    n_nodes : int
        Number of nodes (e.g. 264 for NFA network).
    edges_or_adj : list of edge tuples (u, v), or boolean adjacency matrix (N, N), or adjacency list.

    Returns
    -------
    dict with:
        'degree_mean': float
        'degree_median': float
        'degree_std': float
        'degree_max': int
        'degree_sequence': List[int]
        'components': List[List[int]]
        'num_components': int
        'largest_component_size': int
        'largest_connected_component_fraction': float
        'isolated_node_count': int
        'isolated_node_fraction': float
    """
    if n_nodes <= 0:
        raise ValueError(f"n_nodes must be positive, got {n_nodes}")

    adj_list: List[List[int]] = [[] for _ in range(n_nodes)]

    if isinstance(edges_or_adj, np.ndarray):
        # Adjacency matrix
        for u in range(n_nodes):
            for v in range(u + 1, n_nodes):
                if edges_or_adj[u, v]:
                    adj_list[u].append(v)
                    adj_list[v].append(u)
    elif isinstance(edges_or_adj, list):
        if len(edges_or_adj) > 0 and isinstance(edges_or_adj[0], list):
            # Already adjacency list
            adj_list = edges_or_adj
        else:
            # List of edge pairs (u, v)
            seen_edges: Set[Tuple[int, int]] = set()
            for edge in edges_or_adj:
                u, v = int(edge[0]), int(edge[1])
                if u != v:
                    pair = (min(u, v), max(u, v))
                    if pair not in seen_edges:
                        seen_edges.add(pair)
                        adj_list[u].append(v)
                        adj_list[v].append(u)

    degrees = np.array([len(adj_list[u]) for u in range(n_nodes)], dtype=np.int32)
    components = compute_connected_components_bfs(n_nodes, adj_list)

    lcc_size = len(components[0]) if components else 0
    lcc_fraction = float(lcc_size / n_nodes)

    iso_count = int(np.count_nonzero(degrees == 0))
    iso_fraction = float(iso_count / n_nodes)

    return {
        "degree_mean": float(np.mean(degrees)),
        "degree_median": float(np.median(degrees)),
        "degree_std": float(np.std(degrees)),
        "degree_max": int(np.max(degrees)),
        "degree_sequence": degrees.tolist(),
        "components": components,
        "num_components": len(components),
        "largest_component_size": lcc_size,
        "largest_connected_component_fraction": lcc_fraction,
        "isolated_node_count": iso_count,
        "isolated_node_fraction": iso_fraction,
    }
