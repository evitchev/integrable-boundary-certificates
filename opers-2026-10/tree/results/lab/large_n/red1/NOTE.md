# RED1 (Fable seat, 2026-10-02) -- operator-level reductions of the Sol 3 / Sol 2 opers at the anchor and at N = 0

Seal `SEAL_RED1.md` 26e50dbb... (03:08:10Z), sent to the lead before any computation.  Everything sealed held.
Starting point: the three-term ODEs of item 307.  Exact rescaling, valid for every k
(yt = x^n/n^(3p^2), tau = T/n, Et = E n^(3p^2-4)):

    Sol 3 (R3):  (tau^2 - PX^2/p^2)(tau^2 - pi^2/p^2) phi = yt^(1+M) (tau - (1+M)/2) phi + Et yt phi,     M = 1/k, p^2 = k/(k+1).
    Sol 2 (R2):  tau (tau^2 - pi^2/p^2) phi = yt^(1+M) ((tau - (1+M)/2)^2 - PX^2/(p^2 k)) phi + Et yt phi.

Nothing in (R3), (R2) is singular at n -> infinity or at n -> 0: in the variable x^n both special points are regular.

## (A) Anchor, t -> 3 (k -> infinity)

OPERATOR LEVEL (exact symbolic, `a_op.py` exit 0; tamper fires).
- (R3) tends coefficientwise to, with xi = -yt, vartheta = xi d/dxi, En = Et:

      L_inf:    (vartheta^2 - PX^2)(vartheta^2 - pi^2) phi = xi (vartheta + 1/2 - En) phi.

- The oscillator -psi'' + (y^2 + (lam^2 - 1/4)/y^2) psi = E psi is, after psi = y^(1/2) e^(-y^2/2) chi(xi = y^2),

      Osc(rho): (vartheta^2 - rho^2) chi = xi (vartheta + 1/2 - En) chi,    rho = lam/2 = (l + 1/2)/2,  En = E/4

  (checked: y^2 x Schrodinger = -4 x this form).  So L_inf IS the oscillator equation with the indicial quadratic
  replaced by the product of the two oscillators' quadratics, rho = PX and rho = pi, i.e. l + 1/2 = 2 PX and 2 pi.
- A1: L_inf has terminating solutions exactly at En = N + 1/2 + rho, rho in {+PX, -PX, +pi, -pi}, N = 0..3 (symbolic), and
  none at generic En: its polynomial spectrum is the UNION of the two oscillator spectra.
- A2: its Frobenius solutions are xi^rho 1F3(rho + 1/2 - En; 2rho+1, rho-pi+1, rho+pi+1; xi); the oscillator's are
  xi^rho 1F1(rho + 1/2 - En; 2rho+1; xi): the same upper parameter.  E enters each Frobenius solution of L_inf only
  through the oscillator's parameter.
- What it is NOT: the direct sum (LCLM) of the two oscillator operators (different irregular type at infinity).  In
  Mellin space it is their multiplicative convolution with a Bessel factor: P_X P_pi phihat(nu) = L phihat(nu-1),
  L/(P_X P_pi) = [L/P_X][1/P_pi] = [1/P_X][L/P_pi].
- NOT done at operator level: the global connection problem of L_inf (Meijer-G connection coefficients are products of
  Gamma functions; I state, but did not compute, that the four Q-functions are 1/Gamma(rho + 1/2 - En)).

CHARGE LEVEL (`a_ch.py` exit 0; control fires).
- A3: the oscillator charges from my engine (n = 2, M = 1) are exactly 2^(2K) B_2K(1/2 + l/2), K = 1..5 -- the Stirling
  coefficients of log Gamma(1/2 + rho - En).
- A4 (a LIMIT, exact rationals at k = 10^4 and 10^8): monic F_(2K-1)(PX, pi; k) -> B_2K(1/2 + PX) + B_2K(1/2 + pi):
  spin 1 exact at every k; spins 3, 5, 7, 9 differ by 6.0e-4, 7.5e-4, 1.5e-3, 3.1e-3 at k = 10^4 and by exactly 10^-4 of
  that at k = 10^8.  Without the 1/2 the difference stays at 2.0 .. 13.1 (control).

## (B) N = 0: Sol 3 at t = 5 (k -> -4), Sol 2 at t = -5 (k -> -3/2)

ANSWER to the lead's hypothesis: the logarithm does NOT come from x^n = 1 + n log x at fixed x.  (Argument, not a
machine check: at fixed x the n -> 0 limit of the equation is d_v^2[(v + e) psi] = -g psi, v = log x, in which E only
shifts the origin.)  And in yt = x^n the family is regular at n = 0.  The mechanism is the BLZ one at M -> 0:

OPERATOR LEVEL (exact operator algebra, `b_op.py` exit 0; tamper fires).  Operators are sums zeta^j f_j(D), D = zeta d/dzeta.
- With zeta = yt^gamma (gamma = (1+M)/3 for Sol 3, 1+M for Sol 2) the leading symbol at large zeta is constant-coefficient
  with characteristic polynomial th^4 - th - Et zeta^(-eps) (Sol 3), th^3 - th^2 - Et zeta^(-eps) (Sol 2), eps = n/(k+1) -> 0.
- At the DOUBLE root (Sol 3: 4 th_0^3 = 1, E_0 = -(3/4) th_0; Sol 2: th_0 = 2/3, E_0 = -4/27), after
  phi = exp(-th_0 zeta/gamma) zeta^(..) chi, i.e. tau -> th_0 zeta + q - gamma D:
    * the zeta^top coefficient is the characteristic polynomial;
    * the next coefficient vanishes IDENTICALLY, for every k and every q (this needs the symmetric ordering
      x^(c/2) L x^(c/2); at a simple root it does not vanish -- tamper);
    * the following coefficient is a pure second-order operator:
        Sol 3:  6 gamma^2 th_0^2 [ (D - d_0)^2 - nu^2 ],   nu^2 = (a^2 + b^2 - gamma^2)/(6 gamma^2),
        Sol 2:  gamma^2 [ (D - d_0)^2 - nu^2 ],            nu^2 = (4 b^2 - 6 a^2 - gamma^2)/(6 gamma^2),
      (a^2, b^2 = the squared X and pi exponents of (R3)/(R2)); at N = 0 BOTH give
        nu^2 = 2 (PX^2 + pi^2) - 1/6 = 4 Delta,   Delta = (PX^2 + pi^2)/2 - 1/24 = -(P^2 + Q^2)/2 (item 261).
- Double-scaling limit eps -> 0, zeta = Z/sqrt(eps), Et = E_0 (1 + eps e - (eps/2) log eps):
  (1/eps)[1 - (Et/E_0) zeta^(-eps)] -> log Z - e (sympy limit), the remaining coefficients carry sqrt(eps), eps, and the
  operator tends coefficientwise to the LOG-RADIAL (Kondo; items 271, 276) operator

        ( (D - d_0)^2 - nu^2 ) chi = g Z^2 (log Z - e) chi,     g = 1/(8 gamma^2 th_0) (Sol 3),  4/(27 gamma^2) (Sol 2).

  The fourth- (third-) order operator reduces to this second-order one on the two colliding branches; the other
  branches vary on the scale Z ~ sqrt(eps) and drop out of the coefficientwise limit (a singular perturbation; the
  decoupling is not proved beyond that).  Only PX^2 + pi^2 survives, as item 261 requires.

CHARGE LEVEL (`b_ch.py` exit 0; limits at h = 1e-8, 1e-16, exact rationals).
- B1 Sol 3, k -> -4: operator charges -> certified Sol-3 table at t = 5 (spins 3, 5, 7, 9: 6.7e-9 .. 1.7e-7, then exactly
  1e-8 of that); the table depends on PX^2 + pi^2 only; and operator - BLZ(M = h) at l^2 = 4 Delta (1+M)^2 is O(h) at all
  four spins.  Spin 9 against the polynomial printed in item 271: 3e-17.
- B2 Sol 2, k -> -3/2: the same against the Sol-2 table at t = -5 (spins 5, 7, 9) and BLZ (spins 3..9); item 271: 5e-16.
So "the c = 1 KdV charges come out of the limit": yes, as the M -> 0 limit of the BLZ family with nu^2 = 4 Delta, which
is the momentum map derived at operator level.

## Operator level versus charge level -- summary
| point | operator level (exact) | charge level |
|---|---|---|
| anchor | coefficientwise limit L_inf; = oscillator form with product indicial polynomial; union of spectra (A1); same hypergeometric parameter (A2) | Stirling/Bernoulli (A3); limit k -> infinity of the charges (A4) |
| N = 0 | coefficientwise double-scaling limit = log-radial Kondo operator with nu^2 = 4 Delta, both solutions | limits k -> -4, -3/2 = tables at t = +-5 = BLZ at M -> 0 |
Limits were taken as limits (two step sizes or sympy limits), never by substitution at the singular value.

## Errors and disclosures
- `a_op.py` v0: my CHECK of the oscillator form substituted derivatives of chi(y^2) by pattern, the pattern did not match,
  and the script exited 2 on a correct identity.  Rewritten in the variable xi with a generic function (v0 and its log kept).
- The seal's control wording for (A) ("drop the 1/2") is what was run.
- Not done: Meijer-G connection coefficients (stated, not computed); a proof that the non-colliding branches decouple;
  the Sol-1 points; the exact point (Opus) and numerics (cc) are other seats' parts.

## Exit codes
a_op real 0, tamper 0 (fires); a_ch 0; b_op real 0, tamper 0 (fires); b_ch 0.

## Files
`SEAL_RED1.md` (+ .sha256, .time), `a_op.py` (+ v0), `a_ch.py`, `b_op.py`, `b_ch.py`, engines `mellin_wkb.py`,
`mellin_gen.py`, `wkb3.py` (byte-identical copies), `cft_F.pkl`, logs, `SHA256SUMS`.
