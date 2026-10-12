# SEAL_VIR5d addendum (post-hoc conjecture, sealed before its hold-out) -- Fable seat, 2026-10-01

After z2_compare.py (real): H-A FAILS as sealed (the antisymmetric power-deformation derivative A has the
right top, and at t = 7/5 the right loss-1 layer, but wrong deeper layers); the full kernel is
two-dimensional only at t = 5/3 (K = 3) and t = 7/5 (K = 5), one-dimensional at 11/5 (K = 4), 17/7 (K = 5)
(sealed prediction (i) refuted); the kernel IS swap-invariant (prediction (ii) holds).

READ OFF the two certified antisymmetric classes (this is the training data; not a test):
    t = 5/3 (i = 1):  (X-Y)[(X-Y)^2 - (X+Y) + 1/4]                    = f(u) f(v),  f(z) = z (z^2 - 1/2)
    t = 7/5 (i = 2):  top (X-Y)^5, loss 1 -(10/3)(X-Y)^3 (X+Y), X^3 coefficient 11/3, X^2 Y coefficient -35/9
                                                                       consistent with f(z) = z (z^2 - 1/3)(z^2 - 4/3)
with u = PX + pi, v = PX - pi, X = PX^2, Y = pi^2.

CONJECTURE P-W.  At the boundary drop points K = 2i+1, k = 1/i, t_i = (2i+3)/(2i+1) (p^2 = 1/(i+1),
c_U = c_W = 13 - 6(i+1) - 6/(i+1) = c_(1,i+1) in the U, W Virasoro pair of VIR2), the kernel of ad I_3 at
weight 4i+2 is two-dimensional, and its PX <-> pi antisymmetric class has vacuum eigenvalue
        prod_{j=-i}^{i} (PX + pi - j p)(PX - pi - j p),        p^2 = 1/(i+1)
(a product of two odd polynomials of degree 2i+1 with equally spaced roots, spacing p = the Virasoro
screening charge beta_+).  In operator variables: prod_j (l0 + l1 - j n)(l0 - l1 - j n) up to scale, i.e. it
vanishes exactly when two of the four dynamical exponents differ by j n, |j| <= i.
TESTS (exact).  T-a: the full certified antisymmetric class at t = 7/5 equals the product (all 15
coefficients; 4 were used above) -- 90%.  T-b (hold-out, nothing seen): i = 3, t = 9/7, K = 7, weight 14:
kernel dimension 2 and antisymmetric class = the product with p^2 = 1/4 -- 75%, if my weight-14 pipeline
is affordable (stated either way).  T-c: at the non-boundary drop points the kernel has NO antisymmetric
class (dimension 1; already seen at 11/5, 17/7; to be checked at 13/7 and 23/9, K = 6, weight 12) -- 85%.
No algebra is named: the reading "product of two singlet currents of weight 2i+1" is a guess and is not
tested here.
