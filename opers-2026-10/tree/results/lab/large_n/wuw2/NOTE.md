# WUW2 -- i = 2 at t = 7/5: the antisymmetric spin-9 class IS W_U^(5) W_W^(5) (and the eigenvalue match is not automatic here)

cc (Opus 5.5 seat), 2026-10-01. Lead's commission. Head only.
- Seal: `SEAL_WUW2.md` 6ef71c92921ae360b16d0dad47d75f28ba4501f12ac522373a91204c0c877b26 (18:02:06). The hash was sent after the command returned.
- No Fable code was opened.
- Exit codes: `wuw2.py` 0 (v1 exit 1, kept in `v1_reduce_bug/`); `k10_full.py` 0.

## Setting
- t = 7/5, p^2 = 1/3, N = -24. Bosons V = U/s, Z = W/s, with s^2 = 2/3, U, W = (phi ± X)/sqrt2, and d d ~ kappa = 3/2. All exponents are rational:
  - pair A e^(2V), e^(2Z); pair B e^(-(2/3)V), e^(-(2/3)Z); Virasoro e^((V+Z)/3), e^(-2(V+Z)).
- T_U = (1/3)(dV)^2 + (2/3) d^2V: c = -7 = c_(1,3) (checked).
- T_total = T_U + T_W - (3/2)(d^2V + d^2Z): c = -23 (checked; = 1 + N).
- Engine: WUW1's code with kappa per contraction (`kwick.py`). Control (0): at kappa = 1 it reproduces WUW1 exactly (I_3; W_U; the unique bilinear (1, 3/2, 1/4)^2 in K_6).

## Construction
- **W^(5) := Res e^(2V)(z) e^(-2V)(w) = S_5(2 d^kV/k!) = (4/15)(dV)^5 + (4/3)(dV)^3 d^2V + (2/3)(dV)^2 d^3V + dV (d^2V)^2 + (1/6) dV d^4V + (1/3) d^2V d^3V + (1/60) d^5V**.
- Verified:
  - T_U-primary of weight 5;
  - annihilated at field level by the short screening e^(-(2/3)V);
  - not in span{d^3 T_U, :T_U dT_U:}.
- The long screening e^(2V) does NOT annihilate it (5 terms).

## Results (all sealed tests PASS)
| test | result |
|---|---|
| GATE I_3 (weight-4 commutant of {e^((V+Z)/3), e^(2V), e^(2Z)}) | dim 1, eigenvalue ∝ cyl(3,4) at 7/5 |
| **P1** eigen(W^(5)) = (1/15) P_V (P_V^2 - 1)(4 P_V^2 - 1), P_V = (sqrt3/2)(P_X + pi) | PASS: roots j/2, j = -2..2, the E_2 string |
| **P2** bilinears f_i(V) g_j(Z) (7 × 7 = 49) lying in K_10 + d | **exactly 1-dimensional = W^(5) ⊗ W^(5)** |
| **P3** G_10 = W^(5) W^(5) in K_10; eigen(G_10) = -(81/4800) E_2 | PASS |
| **P4a** G_10 primary of weight 10 under T_U + T_W | PASS |
| **P4b** NO representative G_10 + dH (H any weight-9 field) is T_total-primary | PASS (consistent with cyl_t2_primary) |
| full K_10 = ker ad I_3, weight 10 mod d (`k10_full.py`) | **dim 2** (kernel 302 = 300 d-exact + 2); certified vev_sol3_w10 at 7/5 in span; G_10 in it; the Sol 3 charge is not ∝ eigen(G_10) |

## Vacuity (stated in advance; one correction)
- The seal said the one-boson weight-5 eigenvalue shapes span 3 dimensions. **Computed: 2** (the 7 monomials give only a 2-dimensional span inside {a P^5 + b P^3 + c P}).
- So at i = 2 the eigenvalue shape is a 1-parameter family up to scale. E_2's shape is one point of it. P1 is a genuine test, unlike WUW1's weight 3.
- The decisive test is P2 (density level, uniqueness among 49 bilinears), and it passes.

## Reconciliation with cyl_t2_primary (the record: 0 primary classes at 7/5 under T_total, T_X, T_Y)
- **No conflict.** G_10 is primary for the Coulomb-gas frame T_U + T_W (c = -7 + -7), the frame in which the screenings of VIR2's Vir_U × Vir_W are weight 1.
- The record's T_total differs by -(3/2)(d^2V + d^2Z). Under T_total no representative of G_10's class is primary (P4b, computed over all H).
- The record tested T_total, T_X and T_Y (the X/Y frame). It did not test the U/W frame.
- So "the extra class is not a primary composite field" holds in those frames and fails in the U/W Coulomb-gas frame, where it is the product of two weight-5 primaries.
- Not done: the record's T_X/T_Y tests in my engine. Its statement is taken as recorded.

## Controls (all pass)
- Antisymmetry of eigen(G_10) under P_X <-> pi.
- Symmetric W_V^2 + W_Z^2: symmetric eigenvalue, NOT in K_10.
- Shape-changing tampers:
  - W + (dV)^2 d^3V/10 leaves K_10 and breaks P1 and P3;
  - S_5 with exponent 21/10 breaks P1 and leaves K_10;
  - E_2 with p^2 + 1/1000 breaks P3.

## Disclosures
- **v1 bug** (`v1_reduce_bug/`): the sparse echelon `reduce` eliminated only the LEADING key. That is correct for span membership (P3 and P4b used only membership) but not a linear normal form, so the relation solver returned dim 0 for P2.
  - The kappa = 1 engine control caught it (it failed, exit 1).
  - Fixed to full reduction; all results re-run. The v1 P3/P4b verdicts were unchanged.
- Read before the seal: the record's cyl_t2_primary notes and the docstring of code/cyl_t2_primary.py.

## Files
- Code: `kwick.py`, `wick.py`, `wuw2.py`, `k10_full.py`.
- Logs: `run_wuw2.log/stdout`, `run_k10_full.log`. JSON: `wuw2.json`.
- Inputs: `inputs/` (vev_sol3_w10, hashed); cyl(3,4) via ~/opus5-w8-clone/code.
