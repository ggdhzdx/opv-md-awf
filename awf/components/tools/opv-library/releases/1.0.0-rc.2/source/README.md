# `tool.opv-library`

Portable helper for source-preserving OPV fragment and molecule
libraries. It reads the library's TSV catalogs, checks coordinate/topology
references and package completeness, searches entries, and emits a reviewable
ztop command plan.

The tool is intentionally conservative:

- ztop fragments must come from one source set because references are
  basename-relative;
- `plan-ztop` never executes ztop;
- checks establish file/reference integrity only, not chemical compatibility,
  force-field quality, morphology, or device performance;
- target-project trajectories, checkpoints and results remain outside AWF.

Run from this source directory with:

```bash
python opv_library.py verify
python opv_library.py search Cb --kind fragment
PYTHONPATH=. python -m unittest discover -s tests -v
```
