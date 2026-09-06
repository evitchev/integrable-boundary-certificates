"""VEV recipe validation (2026-09-02): the cylinder eigenvalue calculator
cyl_vev.py reproduces hep-th/0404195 eqs. (59) and (60) -- the paperclip
<I_3> (6 coefficients) and <I_5> (10 coefficients, with the eq.-(51)
erratum) -- from the transcribed densities of families.py, and the three
cylindrical <I_3> (from mixed_engine.cyl_P4) reduce at N = 1 to the
paperclip (59) at n = -1 (Solutions 1 and 3) and n = 2 (Solution 2),
coefficient by coefficient, as the paper's sec. 3.3 claims for the
densities.  Recipe: read the density on the cylinder (d_z^k X ->
d_w^{k-1} J), cylinder normal ordering on |P,Q> (unpaired J -> momentum,
pair (a,b) -> (-1)^b R^{(a+b)}(0), R = 1/(4 sinh^2(eps/2)) - 1/eps^2),
O(N) cycles -> N, chains -> Q^2, article momenta P_art = i P."""
import sys
from fractions import Fraction as F
import sympy as sp
from cyl_vev import vev_density, vev_coefficients
from families import paperclip_P4, paperclip_P6
from mixed_engine import cyl_P4

fails = []
def require(cond, msg):
    print(("PASS  " if cond else "FAIL  ") + msg)
    if not cond:
        fails.append(msg)

def translate(mono):
    xs = tuple(k for k, s in mono if s == 0)
    ys = [k for k, s in mono if s == 1]
    return (xs, tuple((ys[i], ys[i + 1]) for i in range(0, len(ys), 2)))

n, P, Q = sp.symbols('n P Q')
ART3 = (4+3*n)/(6*(2+n))*P**4 + (2+3*n)/(6*n)*Q**4 + P**2*Q**2 + (3+2*n)/(6*(2+n))*P**2 + (1+2*n)/(6*n)*Q**2 + (11+36*n+18*n**2)/(360*n*(2+n))
ART5 = (-(6+5*n)*(8+5*n)/(120*(2+n)**2)*P**6 - (2+5*n)*(4+5*n)/(120*n**2)*Q**6 - (8+5*n)/(8*(2+n))*P**4*Q**2 - (2+5*n)/(8*n)*P**2*Q**4
        - (4+3*n)*(8+5*n)/(48*(2+n)**2)*P**4 - (2+3*n)*(2+5*n)/(48*n**2)*Q**4 - (2+10*n+5*n**2)/(8*n*(2+n))*P**2*Q**2
        - (30+225*n+250*n**2+76*n**3)/(480*n*(2+n)**2)*P**2 - (28+137*n+206*n**2+76*n**3)/(480*n**2*(2+n))*Q**2
        - (564+3410*n+7385*n**2+5680*n**3+1420*n**4)/(60480*n**2*(2+n)**2))
NORM = {3: lambda nv: nv*(nv+2)/((3*nv+2)*(3*nv+4)), 5: lambda nv: -nv*(nv+2)}

def article_coeffs(expr, W, nv):
    """coefficients of the article polynomial in OUR convention: P_art = iP -> P^{2a} picks (-1)^a"""
    pl = sp.Poly(sp.expand(expr.subs(n, sp.Rational(nv))), P, Q)
    out = {}
    for (i, j), c in zip(pl.monoms(), pl.coeffs()):
        a, b = i // 2, j // 2
        out[(a, b)] = F(str(c)) * (-1) ** (a + b)
    for a in range(W // 2 + 1):
        for b in range(W // 2 + 1 - a):
            out.setdefault((a, b), F(0))
    return out

if __name__ == "__main__":
    # P1: paperclip (59), (60): ALL coefficients (6 and 10) at five n, exactly
    for nv in (F(3), F(7, 3), F(-1, 5), F(5), F(-3, 7)):
        d3 = {translate(m): c for m, c in paperclip_P4(nv).items()}
        d5 = {translate(m): c for m, c in paperclip_P6(nv).items()}
        c3 = {k: v / NORM[3](nv) for k, v in vev_coefficients(d3, F(1), 4).items()}
        c5 = {k: v / NORM[5](nv) for k, v in vev_coefficients(d5, F(1), 6).items()}
        require(c3 == article_coeffs(ART3, 4, nv), f"P1: paperclip <I_3> equals eq. (59) coefficient by coefficient (6 coefficients) at n = {nv}")
        require(c5 == article_coeffs(ART5, 6, nv), f"P1: paperclip <I_5> equals eq. (60) coefficient by coefficient (10 coefficients) at n = {nv}")
    # P2: the paper's erratum matters -- with the printed P_6^{[10]} coefficient (60) is NOT reproduced
    nv = F(3)
    d5 = {translate(m): c for m, c in paperclip_P6(nv).items()}
    printed = dict(d5); key = translate(tuple(sorted([(1, 0)] * 6)))
    printed[key] = printed[key] / nv          # remove the erratum's factor n
    c5p = {k: v / NORM[5](nv) for k, v in vev_coefficients(printed, F(1), 6).items()}
    require(c5p != article_coeffs(ART5, 6, nv),
            "P2: with the printed (uncorrected) eq.-(51) coefficient eq. (60) is NOT reproduced (erratum confirmed independently)")
    # P3: N = 1 anchor -- the cylindrical <I_3> at (N, s) = (1, 0) equals the paperclip (59) at n = -1 (sol 1, 3), n = 2 (sol 2)
    for sol, nn, scale in ((1, F(-1), 6), (3, F(-1), 6), (2, F(2), F(3, 2))):
        dens = cyl_P4(sol, F(1), F(0))
        ours = vev_coefficients(dens, F(1), 4)
        theirs = {k: v * scale for k, v in article_coeffs(ART3, 4, nn).items()}
        require(ours == theirs, f"P3: cylindrical Solution {sol} <I_3> at N = 1 equals {scale} x paperclip eq. (59) at n = {nn}, all 6 coefficients")
    print()
    if fails:
        print(f"VEV VALIDATION: {len(fails)} FAILURE(S)"); sys.exit(1)
    print("CYL VEV VALIDATED (cylinder-normal-ordered momentum-state eigenvalues reproduce hep-th/0404195 eqs. (59)-(60) "
          "coefficient by coefficient, the eq.-(51) erratum is required, and the three cylindrical <I_3> reduce at N = 1 to the paperclip at "
          "n = -1, -1, 2 as the paper states for the densities)")
