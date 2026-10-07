from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from opv_morphology.bulkgraph_cli import analyze_manifest
from opv_morphology.phase_separation import (
    assign_mass_to_grid_cic,
    compute_composition_fields,
    compute_gaussian_filtered_variance,
    compute_structure_factor_scc,
)


class ExecutableAnalysisTests(unittest.TestCase):
    def test_bulkgraph_manifest_runs_contact_and_stacking(self) -> None:
        base = np.array(
            [
                [-0.10, -0.10, 0.0], [0.10, -0.10, 0.0], [0.0, 0.10, 0.0],
                [-0.15, 0.00, 0.0], [0.00, 0.15, 0.0], [0.15, 0.00, 0.0],
                [-0.12, 0.08, 0.0], [0.12, 0.08, 0.0], [0.0, -0.15, 0.0],
            ],
            dtype=np.float64,
        )
        coordinates = np.vstack([base + [2.0, 2.0, 2.0], base + [2.0, 2.0, 2.30]])
        elements = np.array(["C"] * len(coordinates))
        molecules = []
        for molecule_id, start in enumerate((0, 9)):
            molecules.append(
                {
                    "id": molecule_id,
                    "conjugated_indices": list(range(start, start + 9)),
                    "fragments": {
                        "core": list(range(start, start + 3)),
                        "terminal_left": list(range(start + 3, start + 6)),
                        "terminal_right": list(range(start + 6, start + 9)),
                    },
                }
            )
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            np.savez(root / "snapshot.npz", coordinates=coordinates, box=np.array([8.0, 8.0, 8.0]), elements=elements)
            manifest = {
                "schema_version": 1,
                "input_npz": "snapshot.npz",
                "molecules": molecules,
            }
            (root / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
            result = analyze_manifest(root / "manifest.json", overlap_area=False)
        self.assertEqual(result["contact_graph"]["edge_count"], 1)
        self.assertEqual(result["stacking_graph"]["primary"]["edge_count"], 1)

    def test_periodic_cic_conserves_mass_and_zeroes_dc_mode(self) -> None:
        box = np.array([4.0, 4.0, 4.0])
        donor_coords = np.array([[0.1, 0.1, 0.1], [3.9, 3.9, 3.9]])
        acceptor_coords = np.array([[2.0, 2.0, 2.0], [2.2, 2.0, 2.0]])
        donor, _ = assign_mass_to_grid_cic(donor_coords, np.array([12.0, 12.0]), box, (16, 16, 16))
        acceptor, _ = assign_mass_to_grid_cic(acceptor_coords, np.array([14.0, 14.0]), box, (16, 16, 16))
        self.assertAlmostEqual(float(donor.sum()), 24.0, places=10)
        self.assertAlmostEqual(float(acceptor.sum()), 28.0, places=10)
        _, delta, _, fraction = compute_composition_fields(acceptor, donor)
        self.assertAlmostEqual(float(delta.mean()), 0.0, places=12)
        self.assertAlmostEqual(fraction, 28.0 / 52.0, places=12)
        q, scc, _, _ = compute_structure_factor_scc(delta, box, delta.shape)
        self.assertEqual(q.ndim, 1)
        self.assertEqual(scc.shape, q.shape)
        self.assertTrue(np.all(np.isfinite(scc)))
        variances = compute_gaussian_filtered_variance(delta, box, delta.shape, [1.0])
        self.assertGreaterEqual(variances["sigma_1nm"], 0.0)

    def test_bulkgraph_rejects_non_orthorhombic_box(self) -> None:
        from opv_morphology.bulkgraph_core.geometry import pbc_distance_matrix

        with self.assertRaisesRegex(ValueError, "orthorhombic"):
            pbc_distance_matrix(np.zeros((1, 3)), np.ones((1, 3)), np.array([4, 4, 4, 90, 80, 90]))


if __name__ == "__main__":
    unittest.main()
