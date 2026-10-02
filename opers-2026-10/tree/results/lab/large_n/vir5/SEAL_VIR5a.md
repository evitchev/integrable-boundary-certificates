# SEAL_VIR5a -- Fable seat.  The generic-t continuation of the family C_k.
Written before any comparison of the continued operator with certified data at a non-integer k.

## 0. Disclosure
Done before sealing, ODE side only: the Mellin-WKB engine `mellin_wkb.py` and its validation
`e0_validate.py` (exit 0): at INTEGER k = 2, 3, 4 it reproduces, coefficient by coefficient at
spins 1..9, the chain engine of VIR4 (which matched the certified tables), including the
missing spin 7 at k = 3; and BLZ.  A first version had a truncation bug (all derivative terms
lost); found by the validation, fixed, log kept.  NOT done: any run at non-integer k; any
comparison at non-integer k.  The spin-11 and spin-13 Solution-3 tables are open to me since
VIR4d; that is why the whole prediction file is hash-registered before any comparison.

## 1. The operator (no free parameter)
For any real k (k != 0, -4), with n = k + 4, M = 1/k, s = (n-1)/2, theta = x d/dx:
     H(theta) psi = x^n ( x^{n/k} - E ) psi,
     H(theta) = [ (s - theta)^2 - l_0^2 ] [ (s - theta)^2 - l_1^2 ]
                (n/k)^k Gamma( u + (k+1)/2 ) / Gamma( u - (k-1)/2 ),     u = k (s - theta)/n,
     (l_0, l_1) = n sqrt((k+1)/k) (PX, pi),   pi = P_phi - rho.
Fibre: k = p^2/(1 - p^2), p^2 = 2(t-1)/(t+1), i.e. t = (3k+2)/(k+2).  For integer k >= 1 the Gamma
ratio is the product over the frozen exponents and this is C_k (items 294-295).
HOW IT ACTS.  As a Mellin-space difference equation: for psi = int psihat(nu) x^nu d nu,
     H(nu) psihat(nu) = psihat(nu - n - n/k) - E psihat(nu - n).
Equivalently a pseudo-differential operator in v = log x with symbol H.
HOW THE WKB COEFFICIENTS ARE DEFINED.  psi = exp(S), theta S = s - T.  The exact conjugation
identity e^{-S} H(theta) e^{S} . 1 = sum_m Ht^(m)(T) C_m, with
sum_m mu^m C_m = exp( sum_{r>=2} (-1)^(r+1) mu^r T^(r-1)/r! ) and Ht(T) := H(s - T), holds for
polynomial H and is taken as the definition for the Gamma-ratio symbol through its large-T
expansion Ht(T) = T^n (1 - l_0^2/T^2)(1 - l_1^2/T^2) sum_j c_j(k) (n/(kT))^(2j),
c_j(k) = k(k-1)...(k-2j+1) [t^(2j)] ( t/(2 sinh(t/2)) )^(k+1).  With E -> -E and
x = E^(k/n) y:  T = y p^(1/n) (1 + sum_i W_i),  p = 1 + y^{n/k};  the charge polynomial of spin
i - 1 is R_i, defined by  - int_0^infty T_i dy/y = B(alpha, beta_0) (k/n) R_i  (analytic
continuation in the Beta arguments), alpha = -(i-1) k/n, beta_0 = (i-1)(k+1)/n; then
l_i^2 -> n^2 ((k+1)/k)(PX^2, pi^2) and division by the PX^i coefficient ("monic").
RATIONALITY IN k.  Every step is a finite sequence of ring operations with coefficients in
Q(k): n, M, c_j(k) (a polynomial in k of degree 3j), binomial coefficients of n - 2j - m,
the derivative weights a + (n/k) b, and Pochhammer ratios in (i-1)(k+1)/n and (i-1)/n.  So each
monic coefficient is a rational function of k; no interpolation in k is used anywhere.  Degree
bound: I do not derive one a priori; the certified coefficients themselves are rational in
t, hence in k, with denominators (read off the tables) of degree at most 2 x spin, and the
claim tested is equality of the two rational functions at the fibres below.

## 2. Fibres (screen first: table coefficients finite, PX-top nonzero; note spins in n Z)
A. inside the verified range 5/3 < t < 3, non-integer k:
     t = 9/4  (k = 10/3, n = 22/3),   t = 21/11 (k = 5/3, n = 17/3; a hold-out fibre of the record's tables).
B. 1 < t < 5/3  (0 < k < 1), outside the verified range:
     t = 3/2  (k = 2/3, n = 14/3),    t = 7/5 (k = 1/2, n = 9/2; spin 9 = 2n: I predict NO spin-9 charge
     from the operator there, as at k = 1, 3, 5 for spin n).
C. p^2 > 1 (k < -1): t > 3 and t < -1:
     t = 4 (k = -6, n = -2),   t = 11/2 (k = -18/5, n = 2/5),   t = -2 (k = -6/5, n = 14/5).
D. |t| < 1 (N > 25, -1 < k < 0):   t = 1/2 (k = -2/5, n = 18/5).
Control: k = 4 (t = 7/3) from the continued formula must equal the registered VIR4d prediction
at spins 11, 13 (b8446a9f...) and the certified values at spins 1..9.

## 3. Protocol
1. Compute, for every fibre, the monic polynomials of spins 1, 3, ..., 13 from the operator of
   sec. 1 -> PREDICTION_VIR5a.json.  Send its sha256 to the lead.  No comparison before
   "registered".
2. Compare: spin 1 against PX^2 + pi^2 - 1/6; spins 3, 5, 7 against cft_F.pkl (gate-checked);
   spin 9 against vev_sol3_w10.json; spins 11, 13 against the anchor11 / anchor13 tables;
   conversion P_art^2 = -PX^2, Q_art^2 = rho^2 - pi^2.  Count mismatches per (fibre, spin).
PASS(fibre) = every coefficient of every spin agrees exactly, except spins that are integer
multiples of n, where the operator is expected to give no charge.

## 4. Predictions
A: both fibres PASS -- confidence 85%.
B: PASS -- 65%.     C: PASS -- 50%.     D: PASS -- 40%.
If A fails, the generic-t continuation is wrong as stated, and I report the first failing
spin and coefficient.  B, C, D are reported cell by cell and do not set the exit code.

## 5. Controls
k = 4 control (sec. 2).  Tamper: at t = 9/4, the Gamma ratio replaced by u^k (c_j = 0 for
j >= 1) must FAIL; M = 1/k + 1/10 must FAIL.
Exit codes of the comparison: 0 = A passes and controls fire; 2 = A fails; 1 = control failure.
