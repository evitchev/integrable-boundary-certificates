# SUPER3: the sl(2|1) / Suzuki ODE test for Sol 1. Opus 5 (2), 2026-10-01.
Seal: SEAL_SUPER3.md, sha256 0c0b70fa..., sent before any WKB. Head only. Exact Fractions. Inputs, copied from the repo (read-only): INPUTS.SHA256.

## Verdict
1. **Step A = (ii), confirmed by computation, and it is a finding about S1's scope.**
   - Suzuki's printed ODE (quant-ph/0003066 eq. (1)), extended by l(l+1)/x^2 (NOT printed by Suzuki; our extension), rescales to psi'' = [kappa^2 (x^(2M) - 1) + kappa alpha x^(M-1) + lam x^-2] psi.
   - The alpha term is kappa-ODD and linear in alpha = sqrt(X) (up to an affine map).
   - With M = (t+3)/(t-1) (the sealed pencil value, NOT fitted) and the affine dictionary alpha^2 = X + tA, lam = r Y + tL fixed at spin 3, **the TOP layer equals the certified Sol 1 top at spins 3, 5, 7 and 9, at all 12 non-degenerate fibres.** The fibres were 2, 4, 6, 7, -4, 5/2, 9/2, the exact point -1/3 and -1/3 +- 1e-6, and 1000, 1e6.
   - Suzuki's Lambda is NOT the pillow's. In w = x^(2M) it is (w-1) w^(1/M-2), i.e. (eps, n) = (2 - 1/M, 1) = ((t+7)/(t+3), 1), against the pillow's (-4/(t+3), -(t+7)/(t+3)).
   - **So S1's sentence "any scalar oper reproducing the certified top has the pillow's quadratic differentials" (lead's reading, item 252) is FALSE outside S1's derivation scope.** That derivation assumes the momenta enter only at kappa^0, linearly. A kappa^1 term linear in sqrt(X) produces a different Lambda with the same top.
   - Tamper M -> M + 1/7 breaks the top at every fibre (run_super3_tamper.log).
2. **Step B: the lower layers FAIL** (sealed P2, P3; I made no pass prediction).
   - At every fibre the spin-3 constant fails.
   - Spin 5 has 6 of its 9 non-normalised coefficients wrong.
   - Spins 7 and 9 have 10 and 15 mismatches. All are below the top (top mismatches 0).
   - This includes the exact point (M = -2, exact and as a limit) and the paperclip limit (t = 1e3, 1e6; the residuals converge to nonzero values, e.g. the spin-5 Y^2 entry -> 36).
   - So THIS family (printed Suzuki ODE + l-term, affine dictionary) is NOT the Sol 1 oper.
   - The family "Suzuki Lambda + kappa sqrt(X) U + a general kappa^0 V_0" is NOT tested. It is the analogue of the closed pillow + V_0 line, and it is now open: a new scalar inverse problem.
3. **Structure (sealed A5, P5): PASS.**
   - Odd spins: the charge part is alpha-EVEN at all fibres, and the integer-Lambda-power (alpha-odd) part is an EXACT total derivative (checked by exact linear solve).
   - Even spins 2, 4, 6, 8: the charge part is nonzero and alpha-ODD at all fibres; its integer-power part is again an exact total derivative.
   - Compared with VIR2 (w8_spectrum.log, Sol 1: one charge at every weight 2..7, X-odd at odd weight, i.e. spins 2, 4, 6): the presence and parity pattern MATCH. This is a structural comparison only; no certified even-spin tables exist.
## Engine validation (before any comparison was read)
- BLZ positive control (ctrl_blz.log), alpha = 0: the I_3 fit (2 equations, 2 unknowns) plus all 3 I_5 coefficients equal the BLZ vacuum polynomials EXACTLY at M = 5, 7/3, 9/5, 11/3, -2, -5/2.
  - Delta = ((2l+1)^2 - 4M^2)/(16(M+1)) comes out, c = 13 - 6(b^2 + 1/b^2), b^2 = 1/(M+1).
  - M = -1/3 (t = -2) is a Gamma-pole fibre, screened. t = -2 is also skipped in the main run.
- Total derivatives integrate to zero (ctrl_engine.log).
- **BUG FOUND AND FIXED:** v0 of the engine omitted -R_0' in R_1 (the order-kappa^1 Riccati equation is R_0' + 2 R_0 R_1 = alpha U).
  - Symptom: every alpha = 0 charge vanished at lam = 0, and the BLZ control failed.
  - The v0 run (run_super3_v0_R1bug.log/.json) is kept. It also showed 'top matches, lower layers fail'; only the fixed run counts.
  - A second v0 artefact, integer-f classes fed into Beta, was removed: those classes are now checked as exact total derivatives.
## Evidence for the scope claim (grep of the author's research record unless stated; for the lead to verify)
| closure | where | fibres | how the momenta / kappa enter | covers kappa sqrt(X) U at generic t? |
|---|---|---|---|---|
| S1 (item 252) | l.21114-21130; my SECTOR1 NOTE | all t | kappa^2 Lambda + X v_X + Y v_Y + V_0: momenta at kappa^0, linear | NO: the derivation assumes the kappa^0 form |
| CURVE1/1c, SFI1' | memory curve1-exact-nonlinear-inverse; items 2xx | 6 + 1 fibres | kappa^0 classes on the pillow | NO |
| 38(ad) kappa-linear W | l.4007-4009 | Sol 3, t = 13/6, two-layer fit; "hit a Beta pole at t = 11/7 ... to be run at a pole-free t" | kappa-linear, half-integer Lambda powers | NO, and NEVER RUN for Sol 1. Memory subleading-layer-exact: "Solution 1 has n - eps = -1, so kappa-odd terms hit Beta poles for it at every t" (in the pillow frame) |
| 38(n) open list | l.2964 | - | "Open: a kappa-linear term" | open entry |
| 38(ah) add.1 | l.4567 | - | "Still open in the scalar class: kappa-odd terms (D-3)" | open entry |
| 38(ah)(x) Codex | l.4604-4613 | the exact point, first order in h | an ANALYTIC h kappa Q source; odd-spin response a total derivative | NO: "Outside the test: NON-analytic deformations sqrt(h) kappa Q" |
| anchor kappa columns | l.1415-1421 (10:15) | the anchor, first order | c kappa (1+u)^(-1/2) g(u), g momentum-INDEPENDENT | NO |
| 38(eh)/(ei) | l.14291-14362 | the exact point | first-order kappa-odd sources (poles in the kappa-odd sector) | NO (perturbative) |
| 38(el) parity theorem | l.14481-14520 | any, perturbative | kappa-odd source at odd orders gives a total derivative; even orders genuine (Kimi's toy q_1^2) | NO: it is the mechanism by which alpha^2 enters, not a closure |
| 38(di) | l.12585 | M2 ladder | kappa-odd surviving scalar source | NO (matrix ladder, perturbative) |
| SEMI5/5b/5c | item 30-31, l.16186ff; todo l.194 | large N, t = -1 + 12/N, first quantum order | kappa e, e = Lambda^(1/2) sum_j (e_j^X X + e_j^Y Y) w^j/N: momentum-LINEAR in X, Y | NO: linear in X, not sqrt(X), and large-N perturbative |
| FL (138) placement | l.16509 | literature | one kappa-odd term in print | not a closure |
## Literature
- Babenko-Smirnov, arXiv:1706.03349 = Nucl. Phys. B 924 (2017) 406. The number is consistent in two printed bibliographies (Tanabe 2604.14899 [28]; Ito-Zhu 2206.08024 [15]); not fetched.
  - Per those printed passages, it concerns the N = 1 SCFT (SU(2)_2 coset / osp(2|2)^(2)), NOT N = 2.
- Tanabe (4.56)-(4.58): one-parameter N = 1 NS values, which cannot be compared with Sol 1's two momenta.
- Kojima and Dorey-Tateo N = 2: NOT located; the arXiv rule forbids fetching.
- No printed N = 2 local-IM vacuum eigenvalues were located.
## Deviations
- Dictionary at spin 3 (no certified spin-1 table; accepted by the lead).
- At fibres without a certified spin-3 point (5/2, 9/2, -1/3, the -1/3 +- 1e-6 limits, 1e3, 1e6), the dictionary was fixed at spin 5 (XY^2-type top, X^2, XY), so spin 3 is untested there. Disclosed in the log as 'fixed at spin 5'.
- The engine bug above.
