# Fast-path kernel records (lab records, NOT certificates)

Retained outputs of `lab/fast_kernel2.py` (modular linear algebra + CRT
lifting of the canonical reduced class basis + exact-identity
verification of every genuine vector in the certified engine).  Each
log ends in one `sol=... dim=... primes=...` line per point.  Under the
certification doctrine these are *exact fast-path computations pending
independent certification*: the registered pure-Python certificates
(`cyl_spin13_sol3_sporadic`, `cyl_spin15_ladder_probe`,
`cyl_spin15_stratification`, `cyl_spin17_double_prediction`) are the
records of record.  Later runs write JSON records with the exact
witness vectors (`FASTK_RECORD=path`).

| log | points | result |
|---|---|---|
| w14_anchors.log | weight 14: sheet 3, 1 at (-143/4,-189/4) | 2, 1 (matches cyl_spin13_sol3_sporadic) |
| w16_t53_merged_sheets.log | weight 16: sheets 2, 3 at (-25/2,-45/2) | 2, 2 (matches the exact kernels) |
| w16_anchors.log | weight 16: four ladder-probe points | 1 each (matches cyl_spin15_ladder_probe) |
| w18_P1_P2.log | weight 18: (3, t4), (2, t5) | 2, 2 -- the double prediction |
| w18_controls.log | weight 18: nine preregistered controls | 1 each |
