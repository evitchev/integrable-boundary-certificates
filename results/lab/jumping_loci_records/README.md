# Lab records of the jumping-loci paper

The records behind the results that the paper *Jumping loci of the cylindrical commuting families: arithmetic ladders on
the curve s^2 = (N-25)(N-1)* (E. Vitchev, 2026) labels [lab result] and that are not rows of its Table 1 (the paper's
Appendix C lists them). Each subfolder holds one computation: a README (what it shows, which statements of the paper rest on
it, how to re-run it), the scripts, the records, `run.sh` (re-runs the scripts into a new output directory, never over the
records) and `compare.sh` (compares that output with the records and prints `ALL IDENTICAL` or `DIFFERENCES`). Run them from
the repository root, as each README shows.

| subfolder | the paper's statements | cost of `run.sh` (one core) | needs beyond `requirements.txt` |
|---|---|---|---|
| `classical_locus_checks/` | the symbols over Q(t); Theorem 1 checked at weights 4-30 and the machine checks of Appendix A; the lifting rule's form at odd spin <= 61; the ledger of Appendix B, with an open-tool cross-check of its entries | about 35 s, below 100 MB | none |
| `stratification_w6_w8/` | the classification over Q-bar at weights 6 and 8 | about 70 s, about 75 MB | Python >= 3.9 |
| `spin21_null_cells/` | the two null cells at spin 21 | about 70 min, about 60 GiB per kernel run | python-flint, jq |
| `singlet_eigenvalues/` | the binomial identities and the conditions of the eigenvalue match | about 10 s, about 90 MB | none |
| `q8_density_test/` | the test at q = 8 | about 4 min, about 1 GiB | python-flint |

The factor eigenvalues at q = 4 and q = 6 and the test of the 49 products are the records
`opers-2026-10/tree/results/lab/large_n/wuw1/` and `wuw2/`. The classification over Q-bar at weights 10, 12 and 14 is the
registered certificate `cyl_stratification_w1014` (`code/cyl_stratification_w1014.py`, `results/cyl_stratification_w1014/`).

All five are independent checks made for the paper's revision in October 2026; each was re-run from this tree for release
v1.2.8, with `compare.sh` printing `ALL IDENTICAL`.
