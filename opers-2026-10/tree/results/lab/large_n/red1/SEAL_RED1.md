# SEAL_RED1 -- Fable seat, 2026-10-01.  Operator-level reductions of the Sol 3 / Sol 2 opers at the anchor and at N = 0.
Commissioned by the lead (operator's request).  Written before any computation of this stage.  Hand derivations below
are what is being sealed; every one is then machine-checked.  Charge level and operator level are kept apart.

Starting point (item 307, exact): Sol 3  [ (T^2-l0^2)(T^2-l1^2) - x^(c/2) T x^(c/2) - E x^n ] phi = 0, T = s - theta,
n = k+4, M = 1/k, c = n(1+M) = n/p^2, l = n (PX, pi)/p.   Sol 2  [ T(T^2-l1^2) - x^(c/2)(T^2-l0^2)x^(c/2) - E x^n ] phi = 0,
n = 2k+3, l1 = n pi/p, l0^2 = n^2 PX^2/(p^2 k).
EXACT RESCALING (valid for every k; yt = x^n/n^(3 p^2), tau = T/n, Et = E n^(3p^2-4)):
   Sol 3:  (tau^2 - PX^2/p^2)(tau^2 - pi^2/p^2) phi = yt^(1+M) (tau - (1+M)/2) phi + Et yt phi.                      (R3)
In (R3) nothing is singular at n -> infinity or n -> 0: both limits are regular in the variable yt = x^n.

## (A) ANCHOR, t -> 3 (k -> infinity, p^2 -> 1, M -> 0)
A-op (operator level).  (R3) tends coefficientwise to  (tau^2-PX^2)(tau^2-pi^2) phi = yt (tau - 1/2 + Et) phi; with
   tau = -vartheta (vartheta = yt d/dyt, gauge s_0 = 0) and xi = -yt:
        L_inf :   (vartheta^2 - PX^2)(vartheta^2 - pi^2) phi = xi (vartheta + 1/2 - En) phi,   En = Et.
   The oscillator -psi'' + (y^2 + (lam^2-1/4)/y^2) psi = E psi is, after psi = y^(1/2) e^(-y^2/2) chi and xi = y^2,
        Osc(rho):  (vartheta^2 - rho^2) chi = xi (vartheta + 1/2 - En) chi,   rho = lam/2,  En = E/4.
   So L_inf is the oscillator equation with the indicial quadratic replaced by the PRODUCT of the two oscillators'
   indicial quadratics (rho = PX and rho = pi, i.e. l + 1/2 = 2 PX and 2 pi).  It is NOT the direct sum (LCLM) of the two
   oscillators; in Mellin space it is their multiplicative convolution with a Bessel factor:
   P_X P_pi phihat(nu) = L phihat(nu-1) factorises as [L/P_X][1/P_pi] = [1/P_X][L/P_pi].
   Sealed consequences, to be machine-checked:
   A1  L_inf has terminating (polynomial x power) solutions exactly at En = N + 1/2 + rho, rho in {+-PX, +-pi},
       N = 0, 1, 2, ...: the union of the two oscillator spectra.  (symbolic check N = 0..3)                       95%
   A2  The Frobenius solutions are xi^rho 1F3(rho + 1/2 - En; ...; xi) (upper parameter the oscillator's), so each
       connection coefficient of L_inf carries E only through 1/Gamma(rho + 1/2 - En): the four Q-functions of L_inf
       are the Q-functions of the two oscillators.  (series check; no global connection problem solved)          90%
A-ch (charge level).
   A3  Stirling: the asymptotic series of log Gamma(1/2 + rho - En) has even coefficients B_2K(1/2 + rho); the
       oscillator charges from my engine (n = 2, M = 1) are monic B_2K(1/2 + l/2) x 2^2K.                          90%
   A4  LIMIT of the certified/operator charges: monic F_(2K-1)(PX, pi; k) -> B_2K(1/2 + PX) + B_2K(1/2 + pi) as k -> infinity,
       K = 1..5, shown at k = 10^4 and 10^8 with exact rationals (difference O(1/k)).                             85%
   Not claimed: that the limit k -> infinity and the large-E expansion commute as analytic statements; A4 is the test.

## (B) N = 0: Sol 3 at t = 5 (k -> -4) and Sol 2 at t = -5 (k -> -3/2); n -> 0 in both
SEALED ANSWER to the lead's hypothesis ("x^n = 1 + n log x gives the log potential"): NOT in that form.  At fixed x the
n -> 0 limit is degenerate (it gives d_v^2[(v + e) psi] = -g psi, v = log x, which has no spectral problem), and in the
variable yt = x^n the family is REGULAR at n = 0 (R3 at p^2 = 4/3: exponents 3/4 and 1, no logarithm).  The logarithm
comes from a different mechanism, which I seal as the derivation to be checked:
B-op.  At n = 0 the three terms of the leading symbol are homogeneous of the same degree (sigma-hat = M dP + dL = 0),
   as for BLZ at M = 0.  With zeta = yt^gamma (gamma = (1+M)/3 for Sol 3, gamma = 1+M for Sol 2) the equation is, at
   large zeta, a constant-coefficient one with characteristic polynomial
        Sol 3:  th^4 - th - Et zeta^(-eps) = 0,        Sol 2:  th^3 - th^2 - Et zeta^(-eps) = 0,     eps -> 0 as n -> 0.
   At the double root (Sol 3: th_0^3 = 1/4, E_0 = -(3/4) th_0; Sol 2: th_0 = 2/3, E_0 = -4/27) two of the exponential
   branches collide.  Conjugating phi = exp(-th_0 zeta/gamma) zeta^q chi and taking the double-scaling limit
        eps -> 0,   zeta = Z/sqrt(eps),   Et = E_0 (1 + eps e - (eps/2) log eps),
   the operator, divided by its zeta^2 (Sol 3) resp. zeta^1 (Sol 2) coefficient, tends coefficientwise to the
   LOG-RADIAL (Kondo, items 271/276) operator
        ( (Z d/dZ)^2 - nu^2 ) chi = g Z^2 (log Z - e) chi,       nu^2 = 2 (PX^2 + pi^2) - 1/6 = 4 Delta,
   Delta = (PX^2 + pi^2)/2 - 1/24 = -(P^2 + Q^2)/2 in article momenta (item 261).  The two other branches decouple.
   By hand I have the slope: nu^2 = 2(PX^2 + pi^2) + const for BOTH solutions (Sol 3: (a^2+b^2)/(6 gamma^2) with
   a^2 + b^2 = (3/4)(PX^2+pi^2), gamma = 1/4; Sol 2: (2 b^2/3 - a^2)/gamma^2 with b^2 = pi^2/3, a^2 = -2 PX^2/9, gamma = 1/3).
   The constant -1/6 and the vanishing of the zeta^3 (resp. zeta^2) coefficient after a suitable q are PREDICTIONS
   (from 4 Delta), to be derived by exact operator algebra.                                           Sol 3: 70%, Sol 2: 65%
B-ch (charge level).
   B1  lim_(k -> -4) of the monic Sol-3 operator charges (engine at k = -4 + 1e-8, -4 + 1e-16, spins 3..9) equals the
       certified Sol-3 table at t = 5, depends on PX^2 + pi^2 only, and equals the monic limit M -> 0 of the BLZ
       charges (my engine, n = 2) at l^2 = 4 Delta (1+M)^2; spin 9 also against the polynomial printed in item 271. 85%
   B2  the same for Sol 2 at k -> -3/2 against the Sol-2 table at t = -5.                                            80%
   The link "log-radial equation -> c = 1 KdV charges" itself is the record's item 271 (formal periods through spin 9);
   B1 re-derives the charges from the BLZ family as M -> 0 and does not use item 271 except for the spin-9 cross-check.

## Controls
(A): tamper En -> the wrong shift (vartheta + En instead of vartheta + 1/2 - En ... i.e. drop the 1/2): A1 must fail.
(B): tamper th_0 -> a simple root of the characteristic polynomial: the zeta^3 (resp. zeta^2) coefficient must NOT vanish.
Charge limits: the same test at k = 10^4 against B_2K(PX) + B_2K(pi) (argument without the 1/2) must fail.
Exit codes: a_op.py, a_ch.py, b_op.py, b_ch.py: 0 = as sealed and controls fire, 2 = not as sealed, 1 = control misbehaves.
Limits are demonstrations at two step sizes with exact rationals (stated as such); operator-level statements are exact
symbolic identities or coefficientwise limits.  No nsimplify.
