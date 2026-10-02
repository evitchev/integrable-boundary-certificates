# SEAL_VIR5c -- Fable seat, 2026-10-01.  The generic-t operator and the record's x/sinh(x) law.
Follow-up (i) to VIR5a, commissioned by the lead.  Written before any computation of this stage.

## 0. Disclosure
Read before sealing (read-only): ansatz-note items 180, 185, 186, 188; Codex's archive files
CONVENTION.md, derive_amplitudes.py, profile_reader.py, profile_check.py and the header fields of
AMPLITUDES_sol3_loss1.json (p, q, r, the three shift amplitudes A_(1,j)(k,t), leading_shift).
Seen: results/lab/anchor15/ contains vev_profile_sol3_w16.json (directory listing only; the
file has NOT been opened by this seat, ever).  NOT done: any run of the engine beyond spin 13;
any numerical evaluation of formula F1 below (it was derived by hand; its K = 1 case gives
PX^2 + pi^2 - 1/6, the known spin-1 eigenvalue); any comparison with a Ward amplitude.
Notation: K = spin index (spin 2K-1) -- the record's "k"; k = fibre parameter of the operator,
k = 2(t-1)/(3-t), n = k + 4 = 2(5-t)/(3-t).  They are different letters here on purpose.

## 1. The question and my answer before computing
QUESTION (lead): is the (t/(2 sinh(t/2)))^(k+1) that generates the Gamma-ratio coefficients of
the generic-t operator the same thing as the record's "x/sinh x law" (items 180, 185-188)?
SEALED ANSWER: NO, they are different objects, and the record's law does not need the Gamma
ratio at all.  The record's law is the leading large-K behaviour (K^(2m) at loss m) of the
layers; it comes from the WKB derivative terms on the CLASSICAL curve (its strength K^(2m) is
the double-pole enhancement at the saddle psi_1 = 0, sec. 2).  The Gamma-ratio coefficients c_j
enter as momentum-like insertions without that enhancement: relative order K^(-1) or lower at
every loss.  At loss 1 the whole effect of the Gamma ratio is to replace the constant n - 1 by
n - k = 4 (the number of non-frozen exponents) in F1.
POSITIVE CLAIM (the real content): the operator DERIVES the record's lower layers.  Concretely
F1 (sec. 2) is a closed formula for the loss-1 layer at every spin, and it equals Codex's
Ward-derived loss-1 amplitudes on Solution 3 (item 185) identically in K; and the operator's
layers at losses 1-4 equal the Ward amplitudes (items 180, 185) at spins 15, 17, 19, beyond
every table this seat has used.

## 2. Formula F1 (hand derivation; to be machine-checked against the engine first)
Operator: Ht(T) psi = Lambda psi, Ht(T) = H_0(T)(1 + ct_1/T^2 + ...), H_0 = T^n phi(T^-2),
phi(w) = (1 - l_0^2 w)(1 - l_1^2 w), ct_1 = c_1(k)(n/k)^2 = -((k-1)/24) c2, c2 = n^2 (k+1)/k.
WKB with theta S = s - T: order 1 is a total derivative; order 2, modulo total derivatives,
    int T_2 dv = int G(Lambda) Lambda'^2 dv - ct_1 int Lambda/(T^2 H_1) dv,
    G = (1/24) d^2/dLambda^2 (H_2/H_1),   H_j = H_0^(j)(T_0).
For Lambda = y^n (1 + y^(n/k)) the Beta integrals give EXACTLY
    int Lambda^(-g-2) Lambda'^2 dv = -(c2/(g+1)) int Lambda^(-g) dv
(the same constant c2 that normalises l^2 = c2 (PX^2, pi^2)).  With w = T^-2, D = w d/dw,
psi_1 = n phi - 2 D phi, and Lagrange inversion of w = w_0 phi(w)^(2/n) (second form, whose
Jacobian is psi_1/(n phi)), the coefficient of w_0^(K-1/2), times -(2K-1), in units l^2 = c2 x:
    Q_K = [w^K] phi^nu + (nu/24) [w^(K-1)] phi^nu ( 4 - 2 D log psi_1 - 4 D^2 log psi_1 ) + (loss >= 2),
    nu = (2K-1)/n,   phi = (1 - PX^2 w)(1 - pi^2 w),   psi_1 = n phi - 2 D phi.        (F1)
Without the Gamma ratio ("stripped" operator, c_j = 0 for j >= 1) the 4 is n - 1.
Record variables: X = P_art^2 = -PX^2, Y = Q_art^2 = rho^2 - pi^2, rho^2 = 2/(t^2-1); the record's
loss-1 layer is the degree-(K-1) part of Q_K(-X, rho^2 - Y), i.e. F1's second term at (-X, -Y)
plus rho^2 times the pi^2-derivative of the top, divided by the X^K coefficient binom(nu, K).
Record form (item 180/185): sum_(j=0)^2 A_(1,j)(K,t) K_(K-1)(p+j, q; 1), p = q = -nu,
K_d(p,q;r) coefficient of X^a Y^(d-a) = binom(d,a) (p)_a (q)_(d-a) r^(d-a)/(p)_d.

## 3. Tests, in this order
S1 (ODE side only, gate).  F1 against the Mellin-WKB engine (mellin_wkb.py, sha256 60795426...,
    unchanged from VIR5a): loss-1 layers at K = 1..7, fibres k = 10/3, 2/3, -6; and the stripped
    variant against the VIR5a tamper_gamma run at k = 10/3.  If F1 fails here the derivation is
    wrong; I correct it using the engine only, disclose, and append an amendment (hashed)
    BEFORE step S3.
S2 (registration).  Engine predictions of the monic polynomials at spins 15, 17, 19
    (engine order 20; if order 20 is too slow at a fibre, 16 = spin 15 only, stated) at the
    fibres t = 9/4, 3/2, 4, 2 (k = 10/3, 2/3, -6, 2), and spin 15 at t = 21/11, 11/2, -2, 1/2, 7/3
    -> PREDICTION_VIR5c.json; sha256 to the lead.  Its spins <= 13 must reproduce the registered
    PREDICTION_VIR5a.json (0c12a7e9...) at the fibres they share.
S3 (derivation test).  F1 in record variables == sum_j A_(1,j) K_(K-1)(p+j,q;1) with Codex's
    AMPLITUDES_sol3_loss1.json, exactly, K = 2..40, at t = 9/4, 21/11, 3/2, 4, 11/2, -2, 1/2, 2, 7/3
    (points where an amplitude or the normalisation has a pole are listed and skipped).
S4 (beyond the tables).  Operator layers at losses 1..4, spins 15, 17, 19, from S2, against
    the Ward amplitudes (Sol 3 loss 1: item 185 archive; losses 2-4: item 180 archive),
    conversion as in sec. 2 (substitution, then homogeneous components).
S5 (blind table).  Spin-15 predictions of S2 against results/lab/anchor15/vev_profile_sol3_w16.json,
    all nine fibres, opened only after the lead's "registered".
S6 (the question itself).  (a) exact: loss-1 layer of the stripped operator minus the true one
    = (nu (n - 5)/24) [w^(K-1)] phi^nu / binom(nu,K)   [operator variables], which in the shift
    basis is a single amplitude that grows like K^1 while A_(1,j) ~ K^2.  (b) illustration, not
    gated: at one fibre, the shift-basis amplitudes of (stripped - true) at loss 2, K = 4..10,
    relative to the true amplitudes (expected to fall off like 1/K).

## 4. Predictions
S1: F1 == engine, 90% (a slip in a hand derivation is the main risk).
S3: PASS at every non-singular (K, t), 85%.     S4: PASS, 85%.     S5: all nine fibres PASS, 85%.
S6(a): PASS, 95%.
If S3 fails while S1 passes, the operator and the Ward layers part company above spin 13 --
I report the first K and the residual, and S4/S5 then decide which side the tables support.

## 5. Controls
Engine: spins <= 13 equal the registered VIR5a predictions (S2).  Negative: F1 with 4 -> n - 1
(no Gamma ratio) must FAIL S3; F1 with rho^2 -> 0 must FAIL S3.  Scoring tamper (after the
seal-hash guard): A_(1,1) + 1/1000 must FAIL S3; one Ward amplitude at loss 2 + 1/1000 must FAIL S4.
Exit codes: s3 and s4 comparison scripts 0 = pass and controls fire, 2 = mismatch, 1 = a control
did not fire.  No nsimplify, no lambdify, exact rationals throughout.

## 6. What is NOT claimed
No all-depth proof of the x/sinh x law; no closed formula at loss >= 2 (engine values only);
Solutions 1, 2 not addressed; the argument for "relative order 1/K at every loss" is a counting
argument, tested exactly only at loss 1.
