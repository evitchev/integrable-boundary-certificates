# SECOND2 — Seal (hypotheses, class, counts, pass/fail rule) — FIXED BEFORE COMPUTING

Qwen peer seat, 2026-10-01. Lead commission: item 309 follow-up.
Question: does a SINGLE scalar second-order ODE with a kappa-ODD momentum term (the Suzuki class
that gave Sol 1 its oper, item 306: -psi'' + [x^(2M) + alpha x^(M-1) + l(l+1)/x^2] psi = E psi,
alpha^2 ~ X) also reproduce Sol 2 or Sol 3?

## Hypotheses

H1 (Sol 2): There exist rational functions M(t), and dictionary parameters
    alpha^2 = s_X * X + c_X,    (l+1/2)^2 = s_Y * Y + c_Y,
all scales free (s_X, s_Y, c_X, c_Y, M all to be determined), such that the TOP-layer coefficients
(highest-degree in X after monic normalisation) of the Suzuki WKB local IMs match the certified
Sol 2 VEV tables at spins 3, 5, 7, 9 for fibres t = 2, 9/4, 4.

H2 (Sol 3): Same statement for Sol 3, EXCEPT that a single kappa-sqrt(X) term breaks PX <-> pi
symmetry. Only the DOUBLET AVERAGE (item 309 B1: charges = alpha-even part averaged over +-alpha)
can match. The test is whether the doublet-average top layer matches certified Sol 3 top.

H0 (negative): No solution set exists for Sol 2 or Sol 3 top layer; the Suzuki class does not
extend beyond Sol 1.

## Class definition

- ODE: -psi'' + [x^(2M) + alpha x^(M-1) + l(l+1)/x^2] psi = E psi
  (Suzuki quant-ph/0003066 eq(1) with l-term extension; Riccati R' + R^2 = Q, R = sum_n kappa^(1-n) R_n)
- Dictionary: alpha^2 = s_X * X + c_X,  lam = (l+1/2)^2 - 1/4 = s_Y * Y + c_Y
- M(t): rational function of t, free at each fibre (NOT fixed to Sol 1 formula (t+3)/(t-1))
- WKB charge extraction: exact via suzuki_wkb.py (Riccati recurrence + continued Beta), same as
  closed_dict.py / run_super3b.py. No numeric integration.
- Top layer: for spin k = 2j-1, the coefficient of X^k in the polynomial after monic normalisation
  (divide by coeff of X^k), matching the convention in cert() from closed_dict.py.

## Fibres and spin range

- Sol 2: t = 2, 9/4, 4 (three fibres, as lead suggested)
- Sol 3: t = 2, 9/4, 4 (three fibres)
- Spins tested: 3, 5, 7, 9 (k = 2, 3, 4, 5 in Riccati indexing)
- All scales free at each fibre: {M, s_X, s_Y, c_X, c_Y} = 5 unknowns per fibre
- Fitting strategy (following run_super3b.py): fix dictionary on spin 3 (monomials XY, X, Y, 1 — 4 eqs)
  or spin 5 if no spin-3 certified point exists. The 5th unknown M is fixed by requiring the ODE
  structure (one additional constraint). If 5 unknowns cannot be fixed from one spin, use two spins.

## Counts

- Total comparisons (positive case): 2 solutions x 3 fibres x 4 spins x (num monomials - 1) coefficients
  At each (sol, t, k): ncoeff = len(certified monomials at weight W=2k) - 1 (monic removed)
  From vev JSONs: w6 has 3 monomials (P^6 Q^0, P^4 Q^2, P^2 Q^4 for top; full table has more),
  w8 has ~5, w10 has ~7. Exact count determined at runtime.
- Approximate total: 2 x 3 x 4 x ~5 = ~120 coefficient comparisons (exact count in log)

## Pass/fail rule

PASS for Sol 2 (or Sol 3): A unique (or finite) solution set {M, s_X, s_Y, c_X, c_Y} exists at ALL
three fibres such that ALL non-normalised top-layer coefficients match the certified tables exactly
(rational equality, no tolerance). The same dictionary form must hold across all three fibres
(i.e., M(t) is a single rational function of t, not fibre-specific constants).

FAIL for Sol 2 (or Sol 3): No solution set exists at any fibre, OR solutions exist at some fibres
but not consistently across all three (no single rational M(t) unifies them).

PARTIAL: Top layer matches at one or two fibres but not all three — report exact fibre-by-fibre status.

## Controls

C1 (Sol 1 plant): At t = 2, the procedure must recover the known Sol 1 dictionary:
    alpha^2 = -4M(M+1) X  (i.e., s_X = -4M(M+1), c_X = 0),
    (l+1/2)^2 = (M-1)^2/4 - (M+1) Y  (i.e., s_Y = -(M+1), c_Y = (M-1)^2/4 - 1/4),
    with M = (t+3)/(t-1) = 5 at t=2.
  All spins 3,5,7,9 must match vev_sol1 tables. This validates the engine + pipeline before
  testing Sols 2 and 3.

C2 (tamper): Same as C1 but with M -> M + 1/10. Must FAIL (at least one coefficient mismatch).

## Engine validation gate (GATE 1)

Before ANY Sol 2/Sol 3 comparison: run C1 (Sol 1 plant at t=2). If C1 fails, STOP — the engine
or pipeline is broken; do not proceed to Sols 2/3. Report the failure.

## Method details

1. For each fibre t and solution sol in {1(control), 2, 3}:
   a. Load certified top-layer coefficients from vev_sol{sol}_w{6,8,10}.json
      (conversion as in cert() of closed_dict.py: parse P^i Q^j keys, substitute t value,
       extract X^(i//2) Y^(j//2) term, divide by X^k coefficient for monic normalisation)
   b. Compute WKB charges R[2k] via riccati(M, nmax=10) from suzuki_wkb.py
   c. Extract non-INTEGER-f charge class: poly in (A, L) where A = alpha^2/4... actually
      A corresponds to alpha and L to lam per closed_dict.py convention
   d. Map A -> s_X*X + c_X, L -> s_Y*Y + c_Y, expand, monic-normalise
   e. Set up equations: mapped(k) - certified(k) = 0 for each monomial
   f. Solve with sympy (resultants/Groebner, ulimit -v 20000000) for unknowns
2. Cross-fibre consistency: check if solutions at t=2, 9/4, 4 lie on a common rational M(t)
3. Report exact solution sets (sympy solve output), not just match/mismatch

## Files

- Engine: suzuki_wkb.py (copied from super3, UNMODIFIED)
- Template: closed_dict.py (reference for cert() conversion)
- Data: vev_sol{1,2,3}_w{6,8,10}.json, vev_w4_points.json (all copied to this stage)
- Output: NOTE.md + SHA256SUMS + full logs (run_second2.log, gate1_plant.log, gate1_tamper.log)

## NOT in scope

- Lower layers (only test if top matches — step 2 of commission)
- Full ODE solution (WKB local IMs only, same as SUPER3b)
- Matching beyond spin 9 (no tables available)
- Sol 1 re-verification beyond the plant control
