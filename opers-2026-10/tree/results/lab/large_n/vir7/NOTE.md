# VIR7 (Fable seat, 2026-10-01) -- Solution 1: the Gamma-symbol form of its oper, zero parameters

Seals: `SEAL_VIR7.md` a7bd6970... (21:50:26Z, class U); `SEAL_VIR7_addendum.md` d3697936... (21:56:22Z, hold-out on S1);
prediction `PREDICTION_VIR7.json` 33af7eb7... (22:02:17Z), each sent to the lead before the computation it governs.
STATUS: T2a and T2b (blind spins 11, 13) passed; the lead's comparison came first.

## Verdict

**1. The sealed class U is EXCLUDED at first quantum order.**  Its seven members (X-even symbols
P = (T^2-l1^2)^(beta_P), L = (T^2-l0^2)(T^2-l1^2)^(beta_L)/T^2, strings of the one length a = 1/M) all reproduce the Sol-1
top and all fail the loss-1 Ward layer (K = 2..20, five fibres).  cc's symbol (item 301) fails the same test (control).

**2. The control I expected to fail is the answer.**  The single-string operator

    S1(t):   (T^2 - l1^2) sigma^a R_a((T - l0)/sigma) psi = x^n (x^(n/a) - E) psi,      T = s - theta,
             a = (t-1)/2,  n = (t+3)/2,  M = 1/a,  sigma = nM = (t+3)/(t-1),
             l1 = n pi/p,   l0 = n PX/(p sqrt(a_Y)),   p^2 = 2(t-1)/(t+1),  a_Y = (t-1)/(t+3),

with NO fitted number, reproduces the certified Solution-1 vacuum eigenvalues:
- T1 (loss-1 layer, formula G1): 0 mismatches against Codex's Sol-1 loss-1 amplitudes, K = 2..20, t = 2, 9/4, 3/2, 4, -2
  (166..208 coefficients per fibre).  The one fitted scale came out C = n^2(M+1)/2 at every fibre, which is l1 = n pi/p.
- T2a (all layers, Mellin-WKB engine generalised to odd powers of 1/T): **0 mismatches in 420 coefficients** -- spin 1;
  losses 1..4 against Codex's amplitudes (items 180/185) at spins 3..9; the full tables vev_sol1_w6/8/10 at spins 5, 7, 9;
  fibres t = 9/4, 7/3, -2, 1/2, 4, 3/2.  Tampers: momentum scale x 11/10 -> 339/420 mismatches; string length a + 1/10 ->
  77/80.
- The operator has X-odd charges at every even spin 2..12 (nonzero, odd in PX), as VIR2's Sol-1 spectrum requires.  No
  tables exist for them; they are in the prediction file.
- nZ rule again: no charge where nu = spin/n is an integer and the classical symbol's power is a polynomial -- spin 7 at
  t = 4, spin 9 at t = 3/2, even spin 8 at t = 7/3 (and spin 5 at t = 2).
My sealed "S fails T1: 90%" was wrong.  S1 is item 293's ODE (Suzuki + centrifugal) in Gamma-symbol form; item 293's
failure below the top was a fixed momentum scale (confirmed by the lead in SUPER3's code and by Opus's SUPER3b, which
now matches Sol 1 with a zero-parameter dictionary in the Schrodinger form).  This note is the independent
Gamma-symbol pipeline.

**3. One class for all three solutions ("class U", hand-derived; machine-checked at the BLZ point).**
The three-term equation  [ P(vartheta) + xi^(1/2) L(vartheta) xi^(1/2) + E xi^b ] phi = 0  (vartheta = xi d/dxi) is
equivalent, by dividing the Mellin transform by a product of Gamma functions, to the two-term Gamma-symbol equation
with symbol P(T) prod_(roots r of L) G_a(T - r), a = b/(1-b), M = 1/a, n = deg P + a deg L (up to rescalings of xi, E).

| solution | P | L | b | a = 1/M | n | order of the three-term ODE |
|---|---|---|---|---|---|---|
| 3 | (T^2-l0^2)(T^2-l1^2) | T | p^2 | k | k+4 | 4 |
| 2 | T (T^2-l1^2) | T^2-l0^2 | p^2 | k | 2k+3 | 3 |
| 1 | T^2-l1^2 | T - l0 | p^2/2 | k/(k+2) = (t-1)/2 | (t+3)/2 | 2 (Suzuki + centrifugal) |

with ONE dictionary: l1 = n pi/p, and l0^2 : l1^2 = PX^2/a_Y : pi^2/a_X (the crossed curve).  So at generic t the Sol-3 oper
is a FOURTH-order and the Sol-2 oper a THIRD-order ordinary differential equation with a three-term potential; the
Gamma-symbol operators of items 298, 301 are their Mellin duals.  Checked by machine only where an independent engine
exists: (T^2 - l^2) G_a(T), n = 2 + a, M = 1/a gives the same monic charges as the Schrodinger equation with
M_S = 1 + 2/a (a = 1/2, 5/8, 2; spins 1..7).  At a = 1/2, spin 5 the Gamma form has no charge (nu = 2) and the
Schrodinger form has one: the nZ gaps are a property of the Gamma form, not of the three-term ODE.  NOT checked: a direct
WKB of the third- and fourth-order three-term equations.

**4. Formula G1** (first quantum order of any Gamma-symbol operator; z = 1/T, nu = (2K-1)/n, c2 = n^2(M+1)):

    Q_2K = [z^2K] phi^nu + (nu c2/24)[z^(2K-2)] phi^nu {(n-1) - nu (n L1 - L2 - L1^2)} + nu [z^(2K-2)] phi^nu b(z),
    L1 = D log phi, L2 = D L1,  b(z) = sum over strings c_1(a) sigma^2/(1 - e z)^2,  c_1(a) = -a(a^2-1)/24.

Validated: it reproduces Codex's loss-1 Ward amplitudes of Sol 3, of Sol 2 (from cc's S2(k)) and of Sol 1 (from S1),
K = 2..20, five fibres each -- the loss-1 Ward layer of all three solutions is now derived from an operator.

## Post-hoc (labelled)
For the even n = 2 base the unique first-order symbol correction reproducing the loss-1 layer is
b/c2 = (1 + 2 a_Y)/24 + ((1 - a_Y)/12)(1+u)/(1-u)^2, u = l0^2/T^2: a string pair of length a_Y with Suzuki's M plus a
constant that no sealed frozen block gives.  Not pursued once S1 passed.

## Errors, deviations, disclosures
- v0 of `h1_scan.py` did not test the control S at all (its X scale was restricted to a finite set that missed the
  crossed-curve value) and printed "S fails".  v1 uses the crossed-curve ratio.  Both scripts and logs kept; member
  verdicts identical.
- v0 of `h2_inverse.py` mis-indexed the series by one (log kept); v1 starts at the spin-1 charge.
- `s1_predict.py` v0 asserted a parity on a zero polynomial and died at t = 7/3, order 9 (the vanishing even-spin charge);
  fixed, log kept.
- The class-U equivalence is a hand derivation; rescalings of xi and E and the sign conventions are not spelled out and
  only the BLZ point is machine-checked.
- S1 is a control promoted to the tested object after the seal's own T1; the addendum seal was written before any engine
  run and before any comparison beyond loss 1.
- Exit codes: g0_validate 0; h1_scan real 3 (no member passes; v1 control line reads "S fails T1: False");
  e0_engine_validate 0; s2_compare real 0, tamper_scale 0 (fires), tamper_a 0 (fires); s3_blind real 0, tamper 0 (fires).

## Not claimed
Even-spin charges are predictions without data.  No statement about excited states.  The three-term ODE forms for
Sols 2, 3 are derived, not directly tested.  Nothing is proved identically in t.

## T2b (blind spins 11, 13): PASS
PREDICTION_VIR7.json was hashed (22:02:17Z) before any comparison of engine output with anything.  The lead registered
it and compared first, from its snapshot, against anchor11/vev_profile_sol1_w12.json and anchor13/vev_profile_sol1_w14.json
(opened by no seat before): 0 mismatches at all six fibres, 28 + 36 coefficients each, 384 in all; its t + 1/100
control fails 27/28 and 35/36.  After the lead's "registered" this seat ran `s3_blind.py`: 0 / 384 (exit 0); table tamper
(+1/1000 on one coefficient per cell) gives 12 mismatches (fires).  Sealed predictions: T2a 55%, T2b 85% -- both passed.

## Files
`SEAL_VIR7.md`, `SEAL_VIR7_addendum.md`, `PREDICTION_VIR7.json` (+ .sha256, .time), `g1_lib.py`, `g0_validate.py`,
`h1_core.py`, `h1_scan.py` (+ v0), `h2_inverse.py`, `mellin_wkb.py` (60795426...), `mellin_gen.py`, `e0_engine_validate.py`,
`s1_predict.py`, `s2_compare.py`, `s3_blind.py`, `s1_pred_t*.json`, logs, `SHA256SUMS`.
