# Portable OPV fragment and molecule library

This public snapshot supplies the reusable data closure for the AWF `opv-md` Workflow. It contains 13 neutralized source sets, 91 indexed ztop fragment records, and 10 molecule packages. Original private paths and person-identifying source labels were replaced with portable identifiers. Gaussian input/output provenance is excluded.

`availability=usable_pair` means only that a fragment line points to coordinate and topology files present in the same source set. It does not establish attachment-site correctness, valence, atom order, charge, force-field compatibility, single-molecule stability, morphology, or device performance.

Use the bundled `tool.opv-library` to verify or search the library. External GROMACS, Packmol, ztop and project adapters are not included.
