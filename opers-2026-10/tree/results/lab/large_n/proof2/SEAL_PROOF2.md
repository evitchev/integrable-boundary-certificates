# PROOF2 seal. Opus 5 (2), 2026-10-02. Before any computation.
Read before sealing (definitions only): item 304 (VIR5d Z0); the file lists of the Codex audit archives (codex_wardsol12_582806a: Sol 3 projection_loss2/3.json, loss 4 in portable_replay_final); VIR8's file list.
Labels: PROVED / VERIFIED / NOT DONE, as in PROOF1.
## (1) Rank-drop vanishing at losses 2 and 3 (and 4 if its law is present), Sol 3
- Locus (item 304): nu = (2K-1)/n = i, i a positive integer, 2i <= K - 1 (Sol 3: n = k+4 = 2(t-5)/(t-3), a_X = a_Y = (t-3)/(t-5), so nu = a_Y (K - 1/2)).
- ROUTE A (record side, sealed):
  - The record's loss-m layer is coef_m(a) = L_m(a, K, t) c_a, where L_m is Codex's Ward-derived law ('projection_over_c_a', exact in Q(a, K, t)) and c_a ∝ binom(nu, a) binom(nu, K - a) (the BL top at x = -nu).
  - Rewrite L_m in (a, K, nu) and factor it. At nu = i the top has zero order [a >= i+1] + [K - a >= i+1], which is >= 1 everywhere since K >= 2i+1.
  - PROVE vanishing for all integers 0 <= a <= K - m, all K and all i with 2i <= K - 1: every denominator factor's pole at nu = i must be dominated, by case analysis on the linear factors, symbolic in (a, K, i).
  - This is a theorem about the RECORD's (Ward) charges at losses 2 and 3. The operator side (Gamma form) at loss >= 2 equals the record only as VERIFIED (VIR5a, spins <= 13). So the operator statement at loss >= 2 is labelled VERIFIED + reduction.
- ROUTE B (operator side, attempted only if A succeeds quickly): a generating-function form G_m in Fable's F1' style, fitted per loss and then proved equal to L_m c_a. NOT promised.
- Priors: Route A proved at loss 2: 60%; loss 3: 50%.
## (2) The formal-WKB statement of item 311(1) (the same cB x B(A0, B0+1) for all P, L, a)
- Plan: a reduction argument via log m(nu) - log m(nu - n) (E-independent), plus VERIFIED extension. Fable's archived engines (wkb3.py, mellin_wkb.py, copied out, read-only) are run at NEW parameter points: non-integer a, random roots of L, deg P / deg L = (2,1), (3,2), (4,1).
- Prior: PROVED 15%; the verified extension 80%.
## (3) H1-H3 for Sol 1: NOT attempted unless (1)-(2) finish early.
Caps (ulimit -v 16 GB); full logs.
