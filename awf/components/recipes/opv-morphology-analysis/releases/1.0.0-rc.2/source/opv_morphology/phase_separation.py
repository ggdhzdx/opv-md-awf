#!/usr/bin/env python3
"""
Phase-Separation and Nanoscale Composition Heterogeneity Analysis Engine.
Implements:
1. Mass-weighted component density fields on periodic grids (trilinear CIC).
2. Local mass fraction phi_acceptor(r) and composition fluctuation delta_phi(r).
3. Concentration-concentration structure factor S_cc(q) via 3D FFT.
4. Characteristic domain size d = 2*pi/q* with strict unresolved/box-scale gate.
5. Gaussian-filtered coarse-grained composition variance sigma_phi^2(sigma) at sigma = 1, 2, 3 nm.
6. Two-phase SHA-256 digest verification on source trajectories.
"""

import sys
import os
import json
import hashlib
import time
import argparse
from typing import Tuple, Dict, Any, List, Optional
import numpy as np

# Standard atomic masses in Daltons (atomic mass units)
STANDARD_ATOMIC_MASSES: Dict[str, float] = {
    "H": 1.008,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "F": 18.9984,
    "S": 32.06,
    "CL": 35.45,
    "BR": 79.904,
}


def get_atomic_mass(element: str, name: str = "") -> float:
    """
    Resolve atomic mass from element symbol or atom name.
    Fails closed with ValueError if mass cannot be resolved (no silent fallback).
    """
    el = element.strip().upper() if element else ""
    if el in STANDARD_ATOMIC_MASSES:
        return STANDARD_ATOMIC_MASSES[el]
    cleaned = "".join([c for c in name.strip() if c.isalpha()]).upper()
    if len(cleaned) >= 2 and cleaned[:2] in STANDARD_ATOMIC_MASSES:
        return STANDARD_ATOMIC_MASSES[cleaned[:2]]
    if len(cleaned) >= 1 and cleaned[:1] in STANDARD_ATOMIC_MASSES:
        return STANDARD_ATOMIC_MASSES[cleaned[:1]]
    raise ValueError(
        f"Unknown atom mass: cannot resolve atomic mass for element='{element}', name='{name}'. "
        f"Known elements: {sorted(STANDARD_ATOMIC_MASSES.keys())}. "
        "Failing closed to prevent silent data corruption (no silent carbon fallback)."
    )


def compute_sha256(filepath: str, chunk_size: int = 65536) -> str:
    """Compute SHA-256 hex digest of a file in streaming chunks."""
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while True:
            chunk = f.read(chunk_size)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def get_box_vectors_and_volume(box: np.ndarray) -> Tuple[np.ndarray, float, np.ndarray]:
    """
    Given a box representation, return:
    - h: 3x3 matrix where rows are lattice vectors a, b, c
    - volume: V = |det(h)|
    - inv_h: inverse of h such that s = r @ inv_h gives fractional coordinates
    """
    if box.shape == (3, 3):
        h = box.astype(np.float64)
    elif box.shape == (6,):
        # [lx, ly, lz, alpha, beta, gamma]
        lx, ly, lz = box[:3]
        alpha, beta, gamma = np.deg2rad(box[3:6])
        cx = lx
        cy = 0.0
        cz = 0.0
        bx = ly * np.cos(gamma)
        by = ly * np.sin(gamma)
        bz = 0.0
        ax = lz * np.cos(beta)
        ay = lz * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
        az = np.sqrt(max(0.0, lz**2 - ax**2 - ay**2))
        h = np.array([[cx, cy, cz], [bx, by, bz], [ax, ay, az]], dtype=np.float64)
    elif box.shape == (3,):
        # Orthorhombic [lx, ly, lz]
        h = np.diag(box).astype(np.float64)
    else:
        raise ValueError(f"Unrecognized box shape: {box.shape}")

    vol = float(np.abs(np.linalg.det(h)))
    if vol <= 0.0:
        raise ValueError(f"Non-positive cell volume: {vol}")
    inv_h = np.linalg.inv(h)
    return h, vol, inv_h


def determine_grid_shape(
    box: np.ndarray,
    grid_spacing_nm: float,
    box_in_angstrom: bool = True,
) -> Tuple[int, int, int]:
    """
    Determine integer grid dimensions (Nx, Ny, Nz) for a simulation cell.
    Parameters:
    - box: simulation box representation (3x3 matrix h, 6-vector, or 3-vector).
    - grid_spacing_nm: nominal grid spacing in nanometers (e.g. 0.20 nm or 0.25 nm).
    - box_in_angstrom: whether box dimensions are in Angstroms.
    Returns:
    - grid_shape: (Nx, Ny, Nz) tuple of integers.
    """
    h, _, _ = get_box_vectors_and_volume(box)
    box_lengths = np.linalg.norm(h, axis=1)
    box_lengths_nm = box_lengths / 10.0 if box_in_angstrom else box_lengths
    nx = int(np.round(box_lengths_nm[0] / grid_spacing_nm))
    ny = int(np.round(box_lengths_nm[1] / grid_spacing_nm))
    nz = int(np.round(box_lengths_nm[2] / grid_spacing_nm))
    return (nx, ny, nz)


def assign_mass_to_grid_cic(
    coords: np.ndarray,
    masses: np.ndarray,
    box: np.ndarray,
    grid_shape: Tuple[int, int, int],
) -> Tuple[np.ndarray, float]:
    """
    Assign particle masses to a 3D regular grid using periodic Cloud-in-Cell (CIC / trilinear).
    Parameters:
    - coords: (N, 3) Cartesian coordinates.
    - masses: (N,) atomic masses.
    - box: box representation (3x3 or 6-vector or 3-vector).
    - grid_shape: (Nx, Ny, Nz).
    Returns:
    - grid_mass: (Nx, Ny, Nz) array of accumulated mass in each voxel.
    - voxel_volume: volume of one voxel cell in nm^3.
    """
    h, vol, inv_h = get_box_vectors_and_volume(box)
    nx, ny, nz = grid_shape
    voxel_vol = vol / (nx * ny * nz)

    if len(coords) == 0:
        return np.zeros(grid_shape, dtype=np.float64), voxel_vol

    # Fractional coordinates: s in [0, 1)^3
    s = np.dot(coords, inv_h) % 1.0

    # Scaled coordinates in grid units
    ux = s[:, 0] * nx
    uy = s[:, 1] * ny
    uz = s[:, 2] * nz

    ix0 = np.floor(ux).astype(np.int64) % nx
    iy0 = np.floor(uy).astype(np.int64) % ny
    iz0 = np.floor(uz).astype(np.int64) % nz

    dx = ux - np.floor(ux)
    dy = uy - np.floor(uy)
    dz = uz - np.floor(uz)

    ix1 = (ix0 + 1) % nx
    iy1 = (iy0 + 1) % ny
    iz1 = (iz0 + 1) % nz

    w0x = 1.0 - dx
    w1x = dx
    w0y = 1.0 - dy
    w1y = dy
    w0z = 1.0 - dz
    w1z = dz

    grid_mass = np.zeros(grid_shape, dtype=np.float64)

    # 8 corner combinations
    np.add.at(grid_mass, (ix0, iy0, iz0), masses * (w0x * w0y * w0z))
    np.add.at(grid_mass, (ix1, iy0, iz0), masses * (w1x * w0y * w0z))
    np.add.at(grid_mass, (ix0, iy1, iz0), masses * (w0x * w1y * w0z))
    np.add.at(grid_mass, (ix0, iy0, iz1), masses * (w0x * w0y * w1z))
    np.add.at(grid_mass, (ix1, iy1, iz0), masses * (w1x * w1y * w0z))
    np.add.at(grid_mass, (ix1, iy0, iz1), masses * (w1x * w0y * w1z))
    np.add.at(grid_mass, (ix0, iy1, iz1), masses * (w0x * w1y * w1z))
    np.add.at(grid_mass, (ix1, iy1, iz1), masses * (w1x * w1y * w1z))

    return grid_mass, voxel_vol


def compute_composition_fields(
    grid_mass_acceptor: np.ndarray,
    grid_mass_donor: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, float, float]:
    """
    Compute local acceptor mass fraction field phi_acceptor(r), fluctuation delta_phi(r),
    spatial mean <phi_acceptor>_spatial, and global mass fraction acceptor_mass_fraction_global.
    Returns:
    - phi_acceptor: (Nx, Ny, Nz)
    - delta_phi: (Nx, Ny, Nz) centered by spatial grid mean (mean(delta_phi) == 0)
    - phi_spatial_mean: float, mean(phi_acceptor)
    - acceptor_mass_fraction_global: float, total_mass_nfa / (total_mass_nfa + total_mass_d18)
    """
    total_mass_nfa = float(np.sum(grid_mass_acceptor))
    total_mass_d18 = float(np.sum(grid_mass_donor))
    total_mass = total_mass_nfa + total_mass_d18

    if total_mass <= 0.0:
        raise ValueError("Total mass across both components is zero.")

    acceptor_mass_fraction_global = total_mass_nfa / total_mass
    grid_total = grid_mass_acceptor + grid_mass_donor

    # In void voxels with zero total density, local phi defaults to acceptor_mass_fraction_global
    phi_acceptor = np.where(grid_total > 0.0, grid_mass_acceptor / np.maximum(grid_total, 1e-15), acceptor_mass_fraction_global)
    
    # Spatial mean across grid voxels
    phi_spatial_mean = float(np.mean(phi_acceptor))
    # Centered fluctuation: strictly guarantees mean(delta_phi) == 0 and FFT k=0 mode is 0
    delta_phi = phi_acceptor - phi_spatial_mean

    return phi_acceptor, delta_phi, phi_spatial_mean, acceptor_mass_fraction_global


def compute_reciprocal_grid(
    box: np.ndarray,
    grid_shape: Tuple[int, int, int],
) -> Tuple[np.ndarray, float, float]:
    """
    Compute 3D array of wavevector magnitudes q(kx, ky, kz) in rad/nm,
    nominal spherical shell width dq_shell, and minimum non-zero reciprocal wavevector q_min.
    Returns:
    - q_mag: (Nx, Ny, Nz) array of |q|
    - dq_shell: nominal spherical shell width dq = 2*pi / V^(1/3)
    - q_min: magnitude of shortest non-zero reciprocal lattice vector min |2*pi*n*inv_h^T|
    """
    h, vol, inv_h = get_box_vectors_and_volume(box)
    nx, ny, nz = grid_shape
    leff = vol ** (1.0 / 3.0)
    dq_shell = 2.0 * np.pi / leff

    # Integer wave indices m in [-N/2, N/2)
    mx = np.fft.fftfreq(nx, d=1.0) * nx
    my = np.fft.fftfreq(ny, d=1.0) * ny
    mz = np.fft.fftfreq(nz, d=1.0) * nz

    mx_grid, my_grid, mz_grid = np.meshgrid(mx, my, mz, indexing="ij")
    m_vecs = np.stack([mx_grid, my_grid, mz_grid], axis=-1)

    # q = 2*pi * m @ (inv_h).T
    # r = s @ h => q . r = 2*pi * m . s
    m_flat = m_vecs.reshape(-1, 3)
    q_flat = 2.0 * np.pi * np.dot(m_flat, inv_h.T)
    q_mag_flat = np.linalg.norm(q_flat, axis=1)
    q_mag = q_mag_flat.reshape(grid_shape)

    non_zero_q = q_mag_flat[q_mag_flat > 1e-6]
    q_min = float(np.min(non_zero_q)) if len(non_zero_q) > 0 else dq_shell

    return q_mag, dq_shell, q_min


def compute_structure_factor_scc(
    delta_phi: np.ndarray,
    box: np.ndarray,
    grid_shape: Tuple[int, int, int],
    dq: Optional[float] = None,
    q_max: Optional[float] = None,
) -> Tuple[np.ndarray, np.ndarray, float, float]:
    """
    Compute spherically-averaged structure factor S_cc(q).
    Returns:
    - q_bins: 1D array of shell center wavevectors (rad/nm)
    - s_cc: 1D array of spherically-averaged S_cc(q)
    - dq_shell: shell width used
    - q_min: minimum non-zero reciprocal wavevector
    """
    h, vol, inv_h = get_box_vectors_and_volume(box)
    nx, ny, nz = grid_shape
    ntot = nx * ny * nz

    q_mag, default_dq, q_min = compute_reciprocal_grid(box, grid_shape)
    dq_shell = dq if dq is not None else default_dq

    dphi_fft = np.fft.fftn(delta_phi)
    psd = (vol / (ntot**2)) * np.abs(dphi_fft) ** 2

    q_flat = q_mag.flatten()
    psd_flat = psd.flatten()

    non_zero = q_flat > 1e-6
    q_active = q_flat[non_zero]
    psd_active = psd_flat[non_zero]

    if q_max is None:
        q_max = float(np.max(q_active))

    # Shell index b for each wavevector: round to nearest integer shell center b * dq_shell
    bin_indices = np.maximum(np.round(q_active / dq_shell).astype(np.int64), 1)

    max_bin = int(np.round(q_max / dq_shell))
    valid_bins = (bin_indices >= 1) & (bin_indices <= max_bin)
    bin_indices = bin_indices[valid_bins]
    psd_active = psd_active[valid_bins]

    counts = np.bincount(bin_indices, minlength=max_bin + 1)
    sums = np.bincount(bin_indices, weights=psd_active, minlength=max_bin + 1)

    occupied = np.where(counts[1:] > 0)[0] + 1
    q_bins = occupied * dq_shell
    s_cc = sums[occupied] / counts[occupied]

    return q_bins, s_cc, dq_shell, q_min


def detect_characteristic_peak(
    q_bins: np.ndarray,
    s_cc: np.ndarray,
    dq_shell: float,
    q_min: float,
) -> Dict[str, Any]:
    """
    Pre-agreed resolution gate for characteristic low-q peak:
    - Resolved interior peak: candidate peak must be in a strictly higher shell than q_min,
      strictly exceed S(q_min), and be decreasing on both sides across truly consecutive shells:
      b_i - b_{i-1} == 1 and b_{i+1} - b_i == 1 with S(q_i) > S(q_{i-1}) and S(q_i) > S(q_{i+1}).
    - Box-scale limited: lowest-q shell is the global maximum and monotonically non-increasing
      towards higher q across consecutive shells: S(q_min) >= max(S(q)) - eps_scale and
      S(q_k) >= S(q_{k+1}) - eps_scale. Annotated as box-scale limited, d = null.
    - No resolved characteristic scale: all other cases without an interior peak. Does NOT infer
      d >= L, and does NOT claim the system is homogeneous.
    """
    if len(q_bins) < 3:
        return {
            "status": "no_resolved_characteristic_scale",
            "reason": "insufficient_spectrum_bins",
            "q_star": None,
            "d_nm": None,
            "s_qmin": float(s_cc[0]) if len(s_cc) > 0 else 0.0,
            "q_min": float(q_min),
            "note": "insufficient spectrum bins to identify interior peak",
        }

    # Explicit integer shell indices of the returned bins
    shell_indices = np.round(q_bins / dq_shell).astype(np.int64)

    # Shell index containing q_min
    b_qmin = max(1, int(np.round(q_min / dq_shell)))
    qmin_matches = np.where(shell_indices == b_qmin)[0]
    if len(qmin_matches) > 0:
        j_qmin = int(qmin_matches[0])
    else:
        # Fallback to closest shell to q_min
        j_qmin = int(np.argmin(np.abs(q_bins - q_min)))
    s_qmin = float(s_cc[j_qmin])

    # Find candidate interior peaks strictly above the shell containing q_min
    candidate_peaks = []
    for i in range(len(s_cc)):
        # Candidate peak must be in a strictly higher shell than q_min shell
        if shell_indices[i] <= b_qmin:
            continue
        # Must have adjacent neighbors to evaluate discrete curvature
        if i == 0 or i == len(s_cc) - 1:
            continue
        # Consecutive shell check: adjacent bins must be truly contiguous shells (no gaps)
        if (shell_indices[i] - shell_indices[i - 1] != 1) or (shell_indices[i + 1] - shell_indices[i] != 1):
            continue
        # Strict interior local maximum: decreasing on both sides
        if s_cc[i] > s_cc[i - 1] and s_cc[i] > s_cc[i + 1]:
            # Must strictly exceed the lowest-q boundary shell
            if s_cc[i] > s_qmin:
                candidate_peaks.append((i, float(q_bins[i]), float(s_cc[i])))

    if candidate_peaks:
        # Select the most prominent resolved low-q interior peak
        best_peak = max(candidate_peaks, key=lambda item: item[2])
        q_star = float(best_peak[1])
        d_nm = float(2.0 * np.pi / q_star)

        return {
            "status": "resolved",
            "reason": "robust_interior_peak_detected",
            "q_star": q_star,
            "d_nm": d_nm,
            "peak_s_cc": float(best_peak[2]),
            "s_qmin": s_qmin,
            "q_min": float(q_min),
            "qmin_shell_idx": int(b_qmin),
            "note": "characteristic nanoscale domain peak resolved",
        }

    # No interior candidate peak: check if lowest-q shell is global max and monotonically non-increasing
    # Numerical tolerance is strictly scaled floating error (no empirical threshold)
    eps_scale = float(max(s_qmin, 1.0) * np.finfo(np.float64).eps * 100.0)
    is_global_max = bool(s_qmin >= np.max(s_cc) - eps_scale)

    is_contiguous_decay = True
    is_non_increasing = True
    for k in range(j_qmin, len(s_cc) - 1):
        if shell_indices[k + 1] - shell_indices[k] != 1:
            is_contiguous_decay = False
            break
        if s_cc[k + 1] - s_cc[k] > eps_scale:
            is_non_increasing = False
            break

    has_decay_contrast = any(float(s_qmin - s_cc[k]) > eps_scale for k in range(j_qmin + 1, len(s_cc)))
    is_box_scale = is_global_max and has_decay_contrast and is_contiguous_decay and is_non_increasing

    if is_box_scale:
        return {
            "status": "box_scale_limited",
            "reason": "qmin_shell_global_maximum_monotonically_non_increasing",
            "q_star": None,
            "d_nm": None,
            "s_qmin": s_qmin,
            "q_min": float(q_min),
            "qmin_shell_idx": int(b_qmin),
            "note": "box_scale_limited: lowest-q shell is global maximum with monotonic decay towards high-q; domain scale limited by simulation box",
        }
    else:
        return {
            "status": "no_resolved_characteristic_scale",
            "reason": "no_interior_local_maximum_above_qmin",
            "q_star": None,
            "d_nm": None,
            "s_qmin": s_qmin,
            "q_min": float(q_min),
            "qmin_shell_idx": int(b_qmin),
            "note": "no resolved characteristic scale within observable q range",
        }


def evaluate_multigrid_robustness(
    peak_primary: Dict[str, Any],
    peak_secondary: Dict[str, Any],
    dq_tol: float,
) -> Dict[str, Any]:
    """
    Preregistered numerical robustness check between primary (N=50) and secondary (N=40) grids.
    Peak position must be stable across grid resolutions (|Delta q*| <= dq_tol).
    """
    p_stat = peak_primary.get("status")
    s_stat = peak_secondary.get("status")

    if p_stat == "resolved" and s_stat == "resolved":
        q_p = float(peak_primary["q_star"])
        q_s = float(peak_secondary["q_star"])
        delta_q = abs(q_p - q_s)
        grid_robust = bool(delta_q <= dq_tol)
        return {
            "grid_robust": grid_robust,
            "primary_status": p_stat,
            "secondary_status": s_stat,
            "delta_q_star": delta_q,
            "dq_tolerance": dq_tol,
            "reason": "peak_stable_across_grids" if grid_robust else "peak_shifts_between_grids",
        }
    elif p_stat == s_stat:
        return {
            "grid_robust": True,
            "primary_status": p_stat,
            "secondary_status": s_stat,
            "delta_q_star": 0.0,
            "dq_tolerance": dq_tol,
            "reason": f"concordant_classification_{p_stat}",
        }
    else:
        return {
            "grid_robust": False,
            "primary_status": p_stat,
            "secondary_status": s_stat,
            "delta_q_star": None,
            "dq_tolerance": dq_tol,
            "reason": f"discordant_grid_status_{p_stat}_vs_{s_stat}",
        }


def align_and_average_scc(
    scc_records: List[Any],
    nominal_dq: Optional[float] = None,
    common_qmin_shell: Optional[int] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Robustly average S_cc(q) curves across frames strictly restricted to the contiguous shell
    range starting from the lowest common shell (no interpolation / no extrapolation).
    Each frame uses its own dq_shell = 2*pi/V^(1/3).
    The aggregated q coordinate is the mean of the actual q_b across frames.
    Non-contiguous shells are never treated as adjacent.
    Stops at the first gap; never selects a higher-q segment that would discard the q_min boundary.
    If fewer than 3 contiguous shells exist, returns the available short spectrum without raising
    an exception so peak detection can cleanly return no_resolved_characteristic_scale.
    Returns:
    - common_q: 1D array of mean wavevector centers (rad/nm) for the common contiguous shells
    - scc_mean: 1D array of aligned mean S_cc(q)
    - scc_std: 1D array of aligned std S_cc(q)
    """
    if not scc_records:
        return np.array([]), np.array([]), np.array([])

    # Resolve nominal dq safely without calling min(diff(empty))
    if nominal_dq is None:
        dq_candidates = [float(rec[2]) for rec in scc_records if len(rec) > 2]
        if dq_candidates:
            nominal_dq = float(np.mean(dq_candidates))
        else:
            diff_candidates = [
                float(np.min(np.diff(rec[0])))
                for rec in scc_records
                if len(rec[0]) >= 2
            ]
            if diff_candidates:
                nominal_dq = float(np.median(diff_candidates))
            else:
                return np.array([]), np.array([]), np.array([])

    shell_maps = []
    q_maps = []
    frame_shell_sets = []

    for rec in scc_records:
        q_arr = np.asarray(rec[0], dtype=np.float64)
        s_arr = np.asarray(rec[1], dtype=np.float64)
        # Any empty frame must fail closed (cannot silently drop frames)
        if len(q_arr) == 0 or len(s_arr) == 0:
            return np.array([]), np.array([]), np.array([])
        dq_f = float(rec[2]) if len(rec) > 2 else nominal_dq
        shells = np.round(q_arr / dq_f).astype(np.int64)
        s_map = {int(b): float(val) for b, val in zip(shells, s_arr)}
        q_map = {int(b): float(val) for b, val in zip(shells, q_arr)}
        shell_maps.append(s_map)
        q_maps.append(q_map)
        frame_shell_sets.append(set(s_map.keys()))

    common_shells_set = set.intersection(*frame_shell_sets)
    if not common_shells_set:
        return np.array([]), np.array([]), np.array([])

    # If common_qmin_shell is specified, it MUST be present in all frames; never fallback to higher shell
    if common_qmin_shell is not None:
        if common_qmin_shell not in common_shells_set:
            return np.array([]), np.array([]), np.array([])
        b_start = int(common_qmin_shell)
    else:
        b_start = min(common_shells_set)

    # Strictly expand contiguous segment from b_start until first gap
    contiguous_shells = [b_start]
    next_b = b_start + 1
    while next_b in common_shells_set:
        contiguous_shells.append(next_b)
        next_b += 1

    s_matrix = np.array(
        [[s_map[b] for b in contiguous_shells] for s_map in shell_maps],
        dtype=np.float64,
    )
    q_matrix = np.array(
        [[q_map[b] for b in contiguous_shells] for q_map in q_maps],
        dtype=np.float64,
    )

    scc_mean = np.mean(s_matrix, axis=0)
    scc_std = np.std(s_matrix, axis=0)
    common_q = np.mean(q_matrix, axis=0)

    return common_q, scc_mean, scc_std


def compute_gaussian_filtered_variance(
    delta_phi: np.ndarray,
    box: np.ndarray,
    grid_shape: Tuple[int, int, int],
    sigmas_nm: List[float] = [1.0, 2.0, 3.0],
) -> Dict[str, float]:
    """
    Compute coarse-grained composition variance sigma_phi^2(sigma) for specified Gaussian filter widths.
    """
    q_mag, _, _ = compute_reciprocal_grid(box, grid_shape)
    dphi_fft = np.fft.fftn(delta_phi)

    variances: Dict[str, float] = {}
    for sigma in sigmas_nm:
        w_sigma = np.exp(-0.5 * (q_mag**2) * (sigma**2))
        filtered_fft = dphi_fft * w_sigma
        filtered_real = np.real(np.fft.ifftn(filtered_fft))
        var_val = float(np.mean(filtered_real**2))
        key = f"sigma_{int(sigma):d}nm" if sigma.is_integer() else f"sigma_{sigma:.1f}nm"
        variances[key] = var_val

    return variances


def analyze_single_frame(
    coords_acceptor: np.ndarray,
    masses_acceptor: np.ndarray,
    coords_donor: np.ndarray,
    masses_donor: np.ndarray,
    box: np.ndarray,
    grid_shape: Tuple[int, int, int],
    sigmas_nm: List[float] = [1.0, 2.0, 3.0],
) -> Dict[str, Any]:
    """
    Analyze one instantaneous frame.
    """
    h, vol, _ = get_box_vectors_and_volume(box)

    grid_nfa, voxel_vol = assign_mass_to_grid_cic(coords_acceptor, masses_acceptor, box, grid_shape)
    grid_d18, _ = assign_mass_to_grid_cic(coords_donor, masses_donor, box, grid_shape)

    phi_acceptor, delta_phi, phi_spatial_mean, acceptor_mass_fraction_global = compute_composition_fields(grid_nfa, grid_d18)

    q_bins, s_cc, dq_shell, q_min = compute_structure_factor_scc(delta_phi, box, grid_shape)
    peak_info = detect_characteristic_peak(q_bins, s_cc, dq_shell, q_min)
    variances = compute_gaussian_filtered_variance(delta_phi, box, grid_shape, sigmas_nm)

    return {
        "volume_nm3": vol,
        "q_min_rad_nm": q_min,
        "phi_spatial_mean": phi_spatial_mean,
        "acceptor_mass_fraction_global": acceptor_mass_fraction_global,
        "q_bins": q_bins.tolist(),
        "s_cc": s_cc.tolist(),
        "dq_shell": dq_shell,
        "peak_analysis": peak_info,
        "filtered_variances": variances,
    }


def run_trajectory_analysis(
    tpr_path: str,
    xtc_path: str,
    out_dir: str,
    start_ps: float = 5000.0,
    end_ps: float = 10000.0,
    cadence_ps: float = 20.0,
    max_frames: Optional[int] = None,
    grid_spacing_primary: float = 0.2,
    grid_spacing_secondary: float = 0.25,
    sigmas_nm: List[float] = [1.0, 2.0, 3.0],
    donor_sel: str = "resname DDN DPI DAC DC1 DC2",
    acceptor_sel: str = "resname ADN AAC ARS ARC",
) -> Dict[str, Any]:
    """
    Streaming trajectory analysis with two-phase SHA256 integrity verification.
    """
    import MDAnalysis as mda

    t0_wall = time.time()
    os.makedirs(out_dir, exist_ok=True)

    # 1. Pre-run SHA-256 verification
    print(f"[Phase-Separation] Verifying pre-run SHA256 of TPR: {tpr_path}")
    tpr_sha256_pre = compute_sha256(tpr_path)
    print(f"[Phase-Separation] Verifying pre-run SHA256 of XTC: {xtc_path}")
    xtc_sha256_pre = compute_sha256(xtc_path)

    # 2. Load universe and resolve selections and masses
    u = mda.Universe(tpr_path, xtc_path)
    ag_donor = u.select_atoms(donor_sel)
    ag_acceptor = u.select_atoms(acceptor_sel)

    if len(ag_donor) == 0:
        raise ValueError(f"donor selection '{donor_sel}' yielded 0 atoms!")
    if len(ag_acceptor) == 0:
        raise ValueError(f"acceptor selection '{acceptor_sel}' yielded 0 atoms!")

    print(f"[Phase-Separation] Selected atoms - donor: {len(ag_donor)}, acceptor: {len(ag_acceptor)}")

    # Extract atomic masses
    masses_donor = np.array([get_atomic_mass(a.element, a.name) for a in ag_donor], dtype=np.float64)
    masses_acceptor = np.array([get_atomic_mass(a.element, a.name) for a in ag_acceptor], dtype=np.float64)

    # 3. Stream frames
    frame_records = []
    frame_scc_primary = []
    frame_scc_secondary = []
    processed_count = 0

    for ts in u.trajectory:
        t_ps = float(ts.time)
        if t_ps < start_ps - 1e-3:
            continue
        if t_ps > end_ps + 1e-3:
            break

        # Check cadence
        dt_rel = (t_ps - start_ps) % cadence_ps
        if dt_rel > 1e-2 and abs(dt_rel - cadence_ps) > 1e-2:
            continue

        box = ts.triclinic_dimensions if ts.triclinic_dimensions is not None else ts.dimensions

        # Grid shapes: determine grid dimensions converting Angstrom box lengths to nm
        grid_shape_p = determine_grid_shape(box, grid_spacing_primary, box_in_angstrom=True)
        grid_shape_s = determine_grid_shape(box, grid_spacing_secondary, box_in_angstrom=True)

        coords_donor = ag_donor.positions / 10.0  # Convert Angstroms to nm
        coords_acceptor = ag_acceptor.positions / 10.0

        box_nm = box.copy()
        if box_nm.shape == (3, 3):
            box_nm = box_nm / 10.0
        elif box_nm.shape == (6,):
            box_nm[:3] = box_nm[:3] / 10.0
        elif box_nm.shape == (3,):
            box_nm = box_nm / 10.0

        # Evaluate on primary grid (0.20 nm)
        frame_res_p = analyze_single_frame(
            coords_acceptor=coords_acceptor,
            masses_acceptor=masses_acceptor,
            coords_donor=coords_donor,
            masses_donor=masses_donor,
            box=box_nm,
            grid_shape=grid_shape_p,
            sigmas_nm=sigmas_nm,
        )
        frame_res_p["time_ps"] = t_ps
        frame_res_p["frame_index"] = ts.frame

        # Evaluate on secondary grid (0.25 nm) for preregistered numerical robustness
        frame_res_s = analyze_single_frame(
            coords_acceptor=coords_acceptor,
            masses_acceptor=masses_acceptor,
            coords_donor=coords_donor,
            masses_donor=masses_donor,
            box=box_nm,
            grid_shape=grid_shape_s,
            sigmas_nm=sigmas_nm,
        )

        frame_records.append({
            "frame_index": ts.frame,
            "time_ps": t_ps,
            "volume_nm3": frame_res_p["volume_nm3"],
            "q_min_rad_nm": frame_res_p["q_min_rad_nm"],
            "phi_spatial_mean": frame_res_p["phi_spatial_mean"],
            "acceptor_mass_fraction_global": frame_res_p["acceptor_mass_fraction_global"],
            "filtered_variances": frame_res_p["filtered_variances"],
            "primary_grid_peak_status": frame_res_p["peak_analysis"]["status"],
            "secondary_grid_peak_status": frame_res_s["peak_analysis"]["status"],
            "q_star": frame_res_p["peak_analysis"]["q_star"],
            "d_nm": frame_res_p["peak_analysis"]["d_nm"],
        })

        frame_scc_primary.append((
            np.array(frame_res_p["q_bins"]),
            np.array(frame_res_p["s_cc"]),
            float(frame_res_p["dq_shell"]),
            float(frame_res_p["q_min_rad_nm"]),
        ))
        frame_scc_secondary.append((
            np.array(frame_res_s["q_bins"]),
            np.array(frame_res_s["s_cc"]),
            float(frame_res_s["dq_shell"]),
            float(frame_res_s["q_min_rad_nm"]),
        ))

        processed_count += 1
        if max_frames is not None and processed_count >= max_frames:
            print(f"[Phase-Separation] Reached max_frames benchmark limit ({max_frames}). Stopping stream.")
            break

    t1_wall = time.time()
    elapsed_sec = t1_wall - t0_wall

    # 4. Post-run SHA-256 verification
    print(f"[Phase-Separation] Verifying post-run SHA256 of TPR: {tpr_path}")
    tpr_sha256_post = compute_sha256(tpr_path)
    print(f"[Phase-Separation] Verifying post-run SHA256 of XTC: {xtc_path}")
    xtc_sha256_post = compute_sha256(xtc_path)

    if tpr_sha256_pre != tpr_sha256_post:
        raise RuntimeError("FATAL: TPR SHA256 changed after execution! Trajectory modified!")
    if xtc_sha256_pre != xtc_sha256_post:
        raise RuntimeError("FATAL: XTC SHA256 changed after execution! Trajectory modified!")

    # 5. Aggregate metrics over processed frames on common contiguous integer shells
    common_q_p, scc_mean_p, scc_std_p = align_and_average_scc(frame_scc_primary)
    common_q_s, scc_mean_s, scc_std_s = align_and_average_scc(frame_scc_secondary)

    # Average variances
    var_keys = [f"sigma_{int(s):d}nm" if s.is_integer() else f"sigma_{s:.1f}nm" for s in sigmas_nm]
    var_summary = {}
    for vk in var_keys:
        vals = [fr["filtered_variances"][vk] for fr in frame_records]
        var_summary[vk] = {
            "mean": float(np.mean(vals)),
            "std": float(np.std(vals)),
            "raw": vals,
        }

    # Time-averaged peak analysis on primary and secondary grids
    # Step size is nominal grid shell dq; base wavevector reference is common_q[0]
    dq_nom_p = float(np.mean([rec[2] for rec in frame_scc_primary])) if frame_scc_primary else 0.5
    dq_nom_s = float(np.mean([rec[2] for rec in frame_scc_secondary])) if frame_scc_secondary else 0.5
    q_min_ref_p = float(common_q_p[0]) if len(common_q_p) > 0 else dq_nom_p
    q_min_ref_s = float(common_q_s[0]) if len(common_q_s) > 0 else dq_nom_s

    # Locate peaks directly on time-averaged spectra using common low-q shell as reference
    avg_peak_primary = detect_characteristic_peak(common_q_p, scc_mean_p, dq_nom_p, q_min_ref_p)
    avg_peak_secondary = detect_characteristic_peak(common_q_s, scc_mean_s, dq_nom_s, q_min_ref_s)

    dq_tol = max(dq_nom_p, dq_nom_s)
    grid_robustness = evaluate_multigrid_robustness(avg_peak_primary, avg_peak_secondary, dq_tol)

    replication_provenance = {
        "status": "provisional_single_replica",
        "replica_robustness_claimed": False,
        "ensemble_resolved": False,
        "note": "Single-replica trajectory observable is provisional. Protocol requires n=3 independent replicas to assess ensemble robustness; no ensemble resolution is claimed from a single replica.",
    }

    # 6. Save outputs
    summary = {
        "execution_context": {
            "tpr_path": tpr_path,
            "xtc_path": xtc_path,
            "tpr_sha256": tpr_sha256_post,
            "xtc_sha256": xtc_sha256_post,
            "processed_frames": processed_count,
            "elapsed_wall_sec": elapsed_sec,
            "sec_per_frame": elapsed_sec / max(1, processed_count),
        },
        "time_averaged_metrics": {
            "phi_spatial_mean": {
                "mean": float(np.mean([fr["phi_spatial_mean"] for fr in frame_records])),
                "std": float(np.std([fr["phi_spatial_mean"] for fr in frame_records])),
            },
            "acceptor_mass_fraction_global": {
                "mean": float(np.mean([fr["acceptor_mass_fraction_global"] for fr in frame_records])),
                "std": float(np.std([fr["acceptor_mass_fraction_global"] for fr in frame_records])),
            },
            "composition_variance": var_summary,
            "peak_analysis_primary_grid": avg_peak_primary,
            "peak_analysis_secondary_grid": avg_peak_secondary,
            "grid_robustness": grid_robustness,
        },
        "replication_provenance": replication_provenance,
        "frame_records": frame_records,
    }

    summary_path = os.path.join(out_dir, "phase_separation_summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    # Save S_cc(q) curve CSV
    csv_path = os.path.join(out_dir, "scc_curve.csv")
    with open(csv_path, "w") as f:
        f.write("q_rad_nm,scc_mean,scc_std\n")
        for q, sm, ss in zip(common_q_p, scc_mean_p, scc_std_p):
            f.write(f"{q:.5f},{sm:.8e},{ss:.8e}\n")

    return summary


def main():
    parser = argparse.ArgumentParser(description="donor/acceptor Phase-Separation Analysis CLI")
    parser.add_argument("--mode", choices=["synthetic_test", "benchmark_run", "full_run"], required=True)
    parser.add_argument("--tpr", type=str, help="Path to input TPR")
    parser.add_argument("--xtc", type=str, help="Path to input XTC")
    parser.add_argument("--out-dir", type=str, default=".", help="Output directory")
    parser.add_argument("--start-ps", type=float, default=5000.0)
    parser.add_argument("--end-ps", type=float, default=10000.0)
    parser.add_argument("--cadence-ps", type=float, default=20.0)
    parser.add_argument("--max-frames", type=int, default=None)
    parser.add_argument("--donor-selection", default="resname DDN DPI DAC DC1 DC2")
    parser.add_argument("--acceptor-selection", default="resname ADN AAC ARS ARC")
    args = parser.parse_args()

    if args.mode == "synthetic_test":
        print("Running synthetic test mode directly via CLI is deprecated; use test_synthetic_phase_separation.py.")
        sys.exit(0)

    if not args.tpr or not args.xtc:
        parser.error("--tpr and --xtc are required for benchmark_run and full_run")

    run_trajectory_analysis(
        tpr_path=args.tpr,
        xtc_path=args.xtc,
        out_dir=args.out_dir,
        start_ps=args.start_ps,
        end_ps=args.end_ps,
        cadence_ps=args.cadence_ps,
        max_frames=args.max_frames,
        donor_sel=args.donor_selection,
        acceptor_sel=args.acceptor_selection,
    )


if __name__ == "__main__":
    main()
