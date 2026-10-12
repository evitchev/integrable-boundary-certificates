# RED2 seal: operator-level relations at the exact point and the paperclip. Opus 5 (2), 2026-10-01. Written before any computation.
Read before sealing: items 299, 303 and 307 of the ansatz note (definitions only); my SUPER3b and SUPER3 NOTEs. Nothing was computed.
Rules: all solver runs capped (ulimit -v 20 GB); full logs; exact arithmetic; every result labelled OPERATOR level or CHARGE level.

## (A) Exact point t = -1/3 (Sol 1) and t = +1/3 (Sol 2)
A1. PENCIL THEOREM (sealed), and a machine check of its invariants.
    - No kappa-independent LOCAL transformation maps Suzuki + centrifugal (any M, alpha != 0) to the pillow. The class is: x = f(zeta) (Moebius, power maps, any f), gauge psi = h(zeta) phi, and Darboux/Crum transformations with kappa-independent seeds.
    - Reason: in Liouville normal form Q -> f'^2 Q(f) + (Schwarzian term), Darboux adds -2 (log phi_0)''. Both preserve the kappa-grading, and the kappa^1 part transforms covariantly: kappa alpha U dx^2 is a quadratic differential, nonzero for alpha != 0. The pillow's kappa-odd part is 0.
    - Also, at M = -2: the local orders of Lambda dz^2 on the 3-point set are {-5/2, -5/2, +1} (Suzuki) against {-1/2, -1, -5/2} (pillow). This rules out every Moebius map even before kappa-parity.
    - Machine: compute both multisets and confirm that the kappa-odd part is covariant under a generic f (symbolic).
A2. NON-LOCAL test.
    - Suzuki at M = -2 times x^4 has polynomial coefficients. Its Laplace transform is a 4th-order ODE, with singular points at p^2 = -E. In s ~ p/sqrt(-E) these become s = +-1, which are the poles of the pillow's quadratic differential written in s = sqrt(w) (pencil).
    - Tests, in this order:
      (i) local exponents of L4 against those of the pillow (gauge s^a (1-s^2)^b allowed);
      (ii) exact right division L4 / (gauged pillow), with unknown gauge exponents and spectral scale;
      (iii) L4 against the symmetric CUBE of the pillow, up to gauge and scale.
    - Prior: a right factor 30%; a symmetric cube 15%; NONE in these classes 55%.
A3. Sol 2 at t = 1/3 (k = -1/2, n = 2, M = -2, c = -2): L3 = T(T^2 - l1^2) - x^(-1)(T^2 - l0^2)x^(-1) - E x^2.
    - Test: is L3 the symmetric SQUARE of a second-order operator? Gauge the T^2 term away, write L3 = d^3 + p d + q, and test q = p'/2 exactly. If yes, Q = -p/4, compared with the pillow and with Suzuki(M = -2).
    - PENCIL PREDICTION: NO. E enters L3 only through the 0th-order term, so it survives only in q while p is E-free. Then q = p'/2 fails, unless an E-dependent rescaling moves E into p; that case is included in the machine test.
## (B) Paperclip t -> infinity
B1. Sol 1 at M -> 1 (pencil, then machine-checked).
    - H_S(alpha) = -d^2 + x^2 + alpha + l(l+1)/x^2 is the radial oscillator at energy E - alpha. Its spectral determinant depends only on a = (2l + 3 - E + alpha)/4, so it is EXACTLY isospectral with the oscillator at E and l' = l + alpha/2.
    - With the SUPER3b dictionary at M = 1 (alpha^2 = -8X, (l+1/2)^2 = -2Y), l' + 1/2 is proportional to sqrt(-X) + sqrt(-Y): ONE member of the rotated c = -2 pair. alpha -> -alpha (Suzuki's paired +-alpha problems) gives the other member.
    - OPERATOR level: the regular solutions are related by an Erdelyi-Kober fractional integral of order alpha/2 in y = x^2 (Kummer M(a, b, z) -> M(a, b + mu, z)). This is a non-local Mellin-convolution transform, not a local one.
    - CHARGE level: the exact log D of ONE member contains alpha-ODD odd-spin pieces. The WKB (and certified) charges are its alpha-EVEN part, i.e. (1/2) sum over the +-alpha doublet. This explains cc's (1/2) sum_+- B_(s+1) formula.
    - Machine: the Stirling expansion of log Gamma(a) against cc's formula, and the alpha-parity split.
B2. Sol 3's k = -2 member (item 299, 4th order) against the rotated pair.
    - Hypothesis: the TENSOR (Hadamard) product L+ (x) L- of the two rotated oscillators (same E; l+- + 1/2 proportional to P_X +- pi), not their LCLM (cc excluded the LCLM).
    - Pencil support:
      - the leading symbol theta^2 (theta^2 - C x(1 - E x)) has branches 0, 0, +-2R, the tensor-product signature;
      - the exponents at 0 of the product in x = y^2 are {1/2 +- P_X, 1/2 +- pi} up to momentum normalisation, against item 299's {1/2 +- l0, 1/2 +- l1};
      - the WKB charges of a tensor product ADD (R = R+ + R-), so it is isospectral with (1/2) sum after monic normalisation.
    - Machine: build L+ (x) L- in theta form, then search gauge x^g, x -> mu x^(+-1), the spectral identification and the momentum scale for exact equality with item 299's operator. Prior: 60% yes.
    - If no: a Mellin test.
## Reporting
Every claim labelled OPERATOR or CHARGE; NOT DONE with a witness where applicable.
