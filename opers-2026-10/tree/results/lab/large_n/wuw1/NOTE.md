# WUW1 -- density-level test of VIR5d's reading at t = 5/3: the antisymmetric class IS W_U W_W

cc (Opus 5.5 seat), 2026-10-01. Lead's commission. Head only.
- Seal: `SEAL_WUW1.md` 425bf68be71c22c620309b09061dd4bea2e0e1909e262d2612e41564d4350fc6 (17:55:46).
- **Procedural error (disclosed):** my first seal message to the lead carried a PLACEHOLDER hash (74e7a7a4..., "17:31:12"). I sent it in the same batch as the hashing command, before the command returned. I corrected it within minutes. The lead registered the correct hash, and the stage then contained only the seal. Nothing was computed before the acknowledgement. A memory rule was added: hash first, send after.
- Read: ~/fable-work/vir5d/NOTE.md only (as allowed). No Fable code was opened.
- Exit codes:
  - `wuw1.py` **1**: one sealed control (TAMPER 2) did not fire, because it is VACUOUS for eigenvalue tests; see below.
  - `diag_tamper2.py` 0 and `posthoc_products.py` 0 (both post-hoc, labelled).

## Setting
- t = 5/3, p^2 = 1/2. Bosons U, W = (phi ± X)/sqrt2.
- In these coordinates the Sol 3 screenings are rational: e^(2U), e^(2W) (pair A); e^(-U), e^(-W) (pair B); e^((U+W)/2) and e^(-2(U+W)) (Virasoro).
- T_U = (1/2)(dU)^2 + (1/2)d^2U has c = -2 (checked: 4th pole -1 = c/2). e^(2U) and e^(-U) have weight 1.

## (1) The weight-3 object (computed; no printed source invoked)
- **W_U = Res_(z->w) e^(2U)(z) e^(-2U)(w) = (4/3)(dU)^3 + 2 dU d^2U + (1/3) d^3U** (the long screening applied to e^(-2U)). Likewise W_W.
- Verified properties:
  - T_U-primary of weight 3;
  - annihilated at FIELD level by the short screening e^(-U);
  - not a T_U-descendant.
- Reported: the long screening e^(2U) maps it to 8(dU)^2 - 8 d^2U, which is nonzero.
- The name "c = -2 Coulomb-gas singlet-type weight-3 current" describes these verified properties. No symplectic-fermion / triplet identification is claimed.

## (2) Results
- **I_3:** the weight-4 commutant of {e^((U+W)/2), e^(2U), e^(2W)} has dim 1, and its cylinder eigenvalue is proportional to the record's cyl(3,4) at t = 5/3 (GATE pass).
- **K_6 = ker ad I_3 at weight 6 mod d: dim 2** (VIR5d's dimension confirmed independently).
  - It contains the weight-6 screening commutant (I_5; dim 1).
  - Its eigenvalues span both certified vev_sol3_w6 and vev_sol2_w6 at t = 5/3 (GATES pass).
- **P2 (sealed): PASS.** G_6 = W_U W_W lies in K_6 (mod d) and is NOT in span{I_5} + d. So it is the second class.
- **P3 (sealed): PASS.** eigen(G_6) = -(2/9) E_1 exactly (from the two printed factorisations: -1/18 vs 1/4).
  - The run gives Sol 2 = (1/6) Sol 3 - (15/4) eigen(G_6) = (1/6) Sol 3 + (5/6) E_1: EXACTLY VIR5d's relation.
- **P1 (sealed): PASS.** eigen(W_U) = (4/3) P_U (P_U^2 - 1/4).
- **The vacuity finding (post-hoc, `diag_tamper2.py`):** the cylinder vacuum eigenvalue of ANY weight-3 field in one boson is const × (P^3 - P/4). Only the (dU)^3 term contributes, and the 1/4 is 3 × (-1/12), the universal cylinder contraction.
  - Hence at i = 1, **P1 and P3 are AUTOMATIC** for any product (U-weight-3) × (W-weight-3). The root spacing 1/2 in P_U is not evidence for the singlet current specifically.
  - (At i >= 2 this is no longer automatic; not tested here.)
  - This is why the sealed TAMPER 2 (exponent 2 -> 21/10, which only rescales the eigenvalue) did not fire. The control failure is recorded: exit 1.
- **The discriminating test (post-hoc, labelled, `posthoc_products.py`):**
  - Among ALL bilinears sum M_ij f_i(U) g_j(W) of weight-3 one-boson fields (f = (dU)^3, dU d^2U, d^3U), those lying in K_6 + d form a ONE-dimensional space, M ∝ (1, 3/2, 1/4) ⊗ (1, 3/2, 1/4).
  - Since W_U = (4/3)(1, 3/2, 1/4), **W_U W_W is the UNIQUE product of weight-3 currents in the kernel**.
  - The exponent-21/10 product and (dU)^3 (dW)^3 are NOT in K_6 + d.
- P4 (reported): the zero mode of G_6 commutes with pair B (e^(-U), e^(-W)) and with both Virasoro screenings, but NOT with pair A (e^(2U)). In this sense the second class is a pair-B-only commutant element.

## (3) Controls
| control | result |
|---|---|
| antisymmetry: eigen(G_6) odd under P_X <-> pi | pass |
| symmetric W_U^2 + W_W^2: eigenvalue symmetric (not antisymmetric) | pass; reported: NOT in K_6 |
| TAMPER 1: G_6 + (dU)^6/100 leaves K_6 and breaks P3 | pass |
| SCORING TAMPER: E_1 with p^2 + 1/1000 | pass (P3 breaks) |
| eigen of total derivatives = 0 (recipe sanity) | pass |
| TAMPER 2: exponent 21/10 | **FAILS as an eigenvalue control (vacuous)**; it fires at density level (post-hoc) |

## Reading
- At t = 5/3 the extra (antisymmetric) spin-5 class is, at DENSITY level, the normal-ordered product W_U W_W of the two c = -2 weight-3 currents of the U and W Coulomb gases. Its eigenvalue is E_1.
- VIR5d's reading holds at i = 1. The eigenvalue match alone is not evidence (automatic at weight 3); the kernel membership with uniqueness is.

## Files
- Code: `wuw1.py`, `diag_tamper2.py`, `posthoc_products.py`, `wick.py` (VIR2-R copy).
- Logs: `run_*.log/stdout`. JSON: `wuw1.json`. Inputs: `inputs/` (vev_sol2_w6, vev_sol3_w6, hashed); cyl(3,4) via ~/opus5-w8-clone/code.
