#!/usr/bin/env python3
"""Small regression suite for the project-local OPV library helper."""

from __future__ import annotations

import contextlib
import csv
import hashlib
import io
import json
import tempfile
import unittest
from argparse import Namespace
from pathlib import Path

import opv_library


class OpvLibraryFixture(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.library = Path(self.tempdir.name) / "library"
        source = self.library / "sources" / "src_test"
        source.mkdir(parents=True)
        (source / "core.top").write_text("topology\n", encoding="utf-8")
        (source / "core.pdb").write_text("coordinates\n", encoding="utf-8")
        (source / "ALL.frg").write_text(
            "D;p=core.top;x=core.pdb;site=D1:1-2;resname=ADN;comment=test core\n",
            encoding="utf-8",
        )
        molecule = self.library / "molecules" / "nfa" / "TestNFA"
        molecule.mkdir(parents=True)
        (molecule / "TestNFA.pdb").write_text("coordinates\n", encoding="utf-8")
        (molecule / "TestNFA.top").write_text("topology\n", encoding="utf-8")
        index = self.library / "index"
        index.mkdir()
        self._write_tsv(
            index / "fragment_catalog.tsv",
            [
                {
                    "fragment_id": "src_test::ALL::line1",
                    "source_id": "src_test",
                    "source_path": "/source/src_test",
                    "category": "nfa",
                    "frg_file": "sources/src_test/ALL.frg",
                    "line": "1",
                    "slot": "D",
                    "top_reference": "core.top",
                    "coordinate_reference": "core.pdb",
                    "resname": "ADN",
                    "site": "D1:1-2",
                    "comment": "test core",
                    "top_present": "True",
                    "coordinate_present": "True",
                    "availability": "usable_pair",
                    "status": "source_present_needs_context_audit",
                }
            ],
        )
        self._write_tsv(
            index / "molecule_catalog.tsv",
            [
                {
                    "molecule_id": "TestNFA",
                    "category": "nfa",
                    "resource_dir": "molecules/nfa/TestNFA",
                    "status": "candidate",
                    "source": "/source/TestNFA",
                }
            ],
        )
        sums = []
        for path in sorted(self.library.rglob("*")):
            if path.is_file() and path.name != "SHA256SUMS":
                sums.append(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(self.library)}")
        (index / "SHA256SUMS").write_text("\n".join(sums) + "\n", encoding="utf-8")

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    @staticmethod
    def _write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    def test_public_default_library_is_discoverable(self) -> None:
        self.assertTrue(opv_library.DEFAULT_LIBRARY.is_dir())
        self.assertEqual(opv_library.DEFAULT_LIBRARY.name, "opv_fragment_library")

    def test_verify_passes_complete_fixture(self) -> None:
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            status = opv_library.command_verify(Namespace(library=self.library))
        self.assertEqual(status, 0)
        report = json.loads(output.getvalue())
        self.assertEqual(report["status"], "PASS")
        self.assertEqual(report["fragment_records_checked"], 1)

    def test_missing_reference_is_a_hard_failure(self) -> None:
        (self.library / "sources" / "src_test" / "core.pdb").unlink()
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            status = opv_library.command_verify(Namespace(library=self.library))
        self.assertEqual(status, 1)
        report = json.loads(output.getvalue())
        self.assertEqual(report["status"], "FAIL")
        self.assertTrue(any("coordinate" in item for item in report["errors"]))

    def test_plan_contains_source_working_directory(self) -> None:
        output = io.StringIO()
        args = Namespace(
            library=self.library,
            fragment_id=["src_test::ALL::line1"],
            build="D",
            output="assembled.top,assembled.gro",
            ztop=Path("/opt/ztop.py"),
            python="python3",
            plan_file=None,
        )
        with contextlib.redirect_stdout(output):
            status = opv_library.command_plan(args)
        self.assertEqual(status, 0)
        plan = json.loads(output.getvalue())
        self.assertEqual(plan["working_directory"], str(self.library / "sources" / "src_test"))
        self.assertEqual(plan["status"], "PLAN_ONLY")

    def test_molecule_requires_coordinate_and_topology(self) -> None:
        row = {"molecule_id": "TestNFA", "resource_dir": "molecules/nfa/TestNFA"}
        self.assertEqual(opv_library.molecule_diagnostics(self.library, row), [])
        (self.library / "molecules" / "nfa" / "TestNFA" / "TestNFA.top").unlink()
        self.assertEqual(
            opv_library.molecule_diagnostics(self.library, row),
            ["no topology file (.top/.itp)"],
        )


if __name__ == "__main__":
    unittest.main()
