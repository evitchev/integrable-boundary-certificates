# High-spin vacancy: final results (2026-08-11, complete)

Route: promoted charges commute with each family's defining low
charge by construction; by Jacobi each remaining pairwise
commutator lies in ker ad_low at an even charge spin; a trivial
kernel at ONE curve point forces a trivial kernel over Q(t) (the
specialization lemma, kernel-checked in
`lean/Specialization.lean`).  The verdicts below certify the
one-point kernels; all 21 ledger instances are therefore PROVED
identically over Q(t).

| sector | verdict | evidence |
|---|---|---|
| cyl sigma=12,14,16 (Sol 1,2,3 each) | VACANT x9, GF(p) tags=d-rank=379/929/2211 | vacancy_cyl_121416_certificate.json (2026-08-24; supersedes the fraction-free vacancy_high batch 2026-08-01/02) |
| A2(2) sigma=18, 20 | VACANT x2, GF(p) tags=d-rank=678/1351 | vacancy_a2_1820_certificate.json (2026-08-25; supersedes the fraction-free/sparse 2026-08-02/03 runs) |
| cyl sigma=18 (Sol 1,2,3) | VACANT x3, GF(p) tags=d-rank=5106 | vacancy_cyl_s18_certificate.{log,json} |
| A2(2) sigma=24 | VACANT, GF(p) tags=d-rank=5117 | vacancy_a2_s24_certificate.{log,json} |
| cyl sigma=20 (Sol 1,2,3) | VACANT x3, GF(p) tags=d-rank=11491/11491/11491 | vacancy_cyl_s20_certificate.{log,json} |

The ledger (21 = 18 + 3):

- Cylindrical, each of Solutions 1, 2, 3 (promoted over the low
  charge I3; no spin-9 gap here): [I5,I7] (s12), [I5,I9] (s14),
  [I5,I11] and [I7,I9] (s16), [I7,I11] (s18), [I9,I11] (s20) --
  18 pair instances.
- A2(2) (promoted over the low charge I5; the family has NO spin-9
  charge -- exponents 6n+-1): [I7,I11] (s18), [I7,I13] (s20),
  [I11,I13] (s24) -- 3 pairs.  Brackets with I5 vanish by
  construction and are not ledger items.

Historical 70--71-point checks stand as independent evidence
(mutual_checks.md); the batch record is
results/vacancy_final_modp.json; the sigma=20 preregistration is
results/vacancy_s20_launch.json.

This document supersedes results/vacancy_high_partial.md.
