#!/usr/bin/env python3
"""Run BulkGraph contact and stacking analysis from a portable NPZ + JSON manifest."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

from .bulkgraph_core.graph import build_contact_graph, build_stacking_graph


def _indices(value: Any, atom_count: int, label: str) -> np.ndarray:
    result = np.asarray(value, dtype=np.int64)
    if result.ndim != 1 or result.size < 3:
        raise ValueError(f"{label} must contain at least three atom indices")
    if np.any(result < 0) or np.any(result >= atom_count):
        raise ValueError(f"{label} contains an out-of-range atom index")
    if len(np.unique(result)) != len(result):
        raise ValueError(f"{label} contains duplicate atom indices")
    return result


def analyze_manifest(manifest_path: Path, frame: int = -1, overlap_area: bool = True) -> dict:
    manifest_path = manifest_path.resolve()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != 1:
        raise ValueError("manifest schema_version must be 1")
    npz_path = (manifest_path.parent / manifest["input_npz"]).resolve()
    if manifest_path.parent not in npz_path.parents:
        raise ValueError("input_npz must stay under the manifest directory")

    with np.load(npz_path, allow_pickle=False) as data:
        coordinates = np.asarray(data[manifest.get("coordinates_key", "coordinates")], dtype=np.float64)
        box = np.asarray(data[manifest.get("box_key", "box")], dtype=np.float64)
        elements = np.asarray(data[manifest.get("elements_key", "elements")]).astype(str)
    if coordinates.ndim == 3:
        coordinates = coordinates[frame]
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError("coordinates must have shape (atoms,3) or (frames,atoms,3)")
    if box.ndim == 2 and box.shape[1] in {3, 6} and box.shape[0] != 3:
        box = box[frame]
    if box.shape not in {(3,), (6,)}:
        raise ValueError("box must be [Lx,Ly,Lz] or [Lx,Ly,Lz,90,90,90]")
    if len(elements) != len(coordinates):
        raise ValueError("elements length must equal atom count")

    molecules = manifest.get("molecules")
    if not isinstance(molecules, list) or len(molecules) < 2:
        raise ValueError("manifest must describe at least two molecules")
    atom_count = len(coordinates)
    conjugated = []
    fragments = []
    used_ids = set()
    for offset, molecule in enumerate(molecules):
        molecule_id = molecule.get("id", offset)
        if molecule_id in used_ids:
            raise ValueError(f"duplicate molecule id: {molecule_id}")
        used_ids.add(molecule_id)
        conjugated_idx = _indices(molecule["conjugated_indices"], atom_count, f"molecule {molecule_id} conjugated_indices")
        conjugated.append(coordinates[conjugated_idx])
        source_fragments = molecule.get("fragments", {})
        mapped = {}
        for public_name, internal_name in (
            ("core", "ADN"),
            ("terminal_left", "AAC_left"),
            ("terminal_right", "AAC_right"),
        ):
            idx = _indices(source_fragments[public_name], atom_count, f"molecule {molecule_id} {public_name}")
            mapped[internal_name] = {
                "coords": coordinates[idx],
                "elements": elements[idx].tolist(),
            }
        fragments.append(mapped)

    protocol = manifest.get("protocol", {})
    contact_cutoff = float(protocol.get("contact_cutoff_nm", 0.45))
    stacking = protocol.get("stacking", {})
    contact = build_contact_graph(conjugated, box, cutoff_nm=contact_cutoff)
    stack = build_stacking_graph(
        fragments,
        box,
        dist_thresh=float(stacking.get("heavy_atom_distance_max_nm", 0.36)),
        angle_thresh=float(stacking.get("tilt_angle_max_deg", 30.0)),
        proj_thresh=float(stacking.get("normal_projection_min_nm", 0.25)),
        compute_overlap_area=overlap_area,
    )
    return {
        "schema_version": 1,
        "input_npz": manifest["input_npz"],
        "frame": frame,
        "molecule_count": len(molecules),
        "contact_graph": contact,
        "stacking_graph": stack,
        "scientific_boundary": (
            "Geometric contacts and stacking metrics describe this supplied frame; "
            "they do not establish equilibrium, electronic coupling, or device performance."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Portable OPV BulkGraph analysis")
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--frame", type=int, default=-1)
    parser.add_argument("--skip-overlap-area", action="store_true")
    args = parser.parse_args()
    if args.output.exists():
        parser.error(f"output already exists: {args.output}")
    result = analyze_manifest(args.manifest, args.frame, not args.skip_overlap_area)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
