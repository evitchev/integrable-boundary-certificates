"""SPEC1 S4 / S2c: exact proof of the E = 0 identity Q^O_j(0) = (m(theta_j)/omega) Q^Gamma_j(0) via Gauss multiplication.  NEW FILE.
Usage: spec1_s2c_gauss.py GAUSS

Inputs (standard facts, cited, not re-derived here):
  (M1) Frobenius coefficients of Meijer G (no resonance, i.e. (H)):  G^{d,0}_{1,d}(z | b; al) = sum_h [prod_{i!=h} Gamma(al_i - al_h)
       / Gamma(b - al_h)] z^{al_h} (1 + O(z)),  G^{D,0}_{0,D}(z | ; ep) = sum_h [prod_{i!=h} Gamma(ep_i - ep_h)] z^{ep_h} (1 + O(z)).
  (M2) Large-z asymptotics (Luke 5.7 / DLMF 16.11): G^{q,0}_{p,q} ~ (2 pi)^((nu-1)/2) nu^(-1/2) exp(-nu z^(1/nu)) z^vt,
       nu = q - p, vt = ((1 - nu)/2 + sum b_j - sum a_j)/nu.
  (GM) Gauss multiplication (DLMF 5.5.6): prod_{i=0}^{k-1} Gamma(u + i/k) = (2 pi)^((k-1)/2) k^(1/2 - k u) Gamma(k u).
These are exactly the formulas spec1_lib.meijer_Q0_ode / meijer_Q0_gamma implement (S2c used them numerically).

Argument.  c = sigma (a + 1), frozen f_i = base + sigma i (i < a), base = s - r + sigma (1 - a)/2, u = (base - theta_j)/c.
  (G1) eps_{f_i} - al_j = u + i/(a+1);  (G2) u + a/(a+1) = beta1' - al_j;  (G3) (a+1) u = g - theta_j/sigma.
  So prod_i Gamma(eps_{f_i} - al_j) = [GM with k = a+1] / Gamma(u + a/(a+1)) = (2 pi)^(a/2) (a+1)^(1/2 - g + theta_j/sigma)
  Gamma(g - theta_j/sigma) / Gamma(beta1' - al_j): the dynamical Gamma products cancel, Gamma(beta1' - al_j) cancels, and
  Q^O_j/Q^G_j = (2 pi)^(-a/2) (a+1)^(g - 1/2 - theta_j/sigma) / Gamma(g - theta_j/sigma) * c^((a+1) al_j) * norm_G/norm_O.
  (G4) sum_i eps_{f_i} = a (s - r)/c;  (G5) X := (d-1) vt_O - D vt_G = 1/2 - g.
  (G6) the log of [Q^O_j/Q^G_j * omega / m(theta_j)] is identically 0 as a linear form in log(2 pi), log sigma, log(a+1),
       log(d-1), log(d+a)  (log c = log sigma + log(a+1)).
Everything is symbolic in s, r, theta, sigma, d and a (a positive integer only enters through GM's k = a + 1 and the finite sum G4).
Numerical side checks: GM itself at 100 digits for k = 2..9 at random u, and the closed form against spec1_lib at two points."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import spec1_lib as SL  # noqa: E402
from spec1_lib import sp, mp  # noqa: E402

FAIL = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""), flush=True)
    if not ok:
        FAIL.append(name)


def symbolic():
    s, r, th, sig, d = sp.symbols("s r theta sigma d")
    a, i = sp.symbols("a i", positive=True, integer=True)
    c = sig * (a + 1)
    base = s - r + sig * (1 - a) / 2
    g = (s - r) / sig + (1 - a) / 2
    al = th / c
    b1 = (s + c / 2 - r) / c
    u = (base - th) / c
    D = d + a
    check("(G1) eps_{f_i} - al_j = u + i/(a+1)", sp.simplify((base + sig * i) / c - al - (u + i / (a + 1))) == 0)
    check("(G2) u + a/(a+1) = beta1' - al_j", sp.simplify(u + a / (a + 1) - (b1 - al)) == 0)
    check("(G3) (a+1) u = g - theta_j/sigma", sp.simplify((a + 1) * u - (g - th / sig)) == 0)
    Sf = sp.summation((base + sig * i) / c, (i, 0, a - 1))
    check("(G4) sum_i eps_{f_i} = a (s - r)/c", sp.simplify(Sf - a * (s - r) / c) == 0, f"{sp.simplify(Sf)}")
    # (G5): (d-1) vt_O = (2-d)/2 + sum al - b1;  D vt_G = (1-D)/2 + sum al + sum eps_f  (sum al over the d dynamical, common)
    SA = sp.Symbol("SA")
    X = ((2 - d) / sp.Integer(2) + SA - b1) - ((1 - D) / sp.Integer(2) + SA + Sf)
    check("(G5) X = (d-1) vt_O - D vt_G = 1/2 - g", sp.simplify(X - (sp.Rational(1, 2) - g)) == 0)
    # (G6) logs as independent symbols; logGamma(g - th/sig) cancels against m and is omitted on both sides
    L2p, Ls, La, Ld1, LD = sp.symbols("L2pi Lsigma La1 Ld1 LD")
    Lc = Ls + La
    logR = (-a / sp.Integer(2) * L2p + ((a + 1) * u - sp.Rational(1, 2)) * La      # Gauss product, inverted
            + (a + 1) * al * Lc                                                     # c^((1-d) al) / c^(-D al)
            + (((D - 1) - (d - 2)) / sp.Integer(2)) * L2p - LD / 2 + Ld1 / 2         # norm_G / norm_O prefactors
            + X * Lc)                                                               # c^(-D vt_G) / c^((1-d) vt_O)
    logm = th / sig * Ls
    logw = (g - sp.Rational(1, 2)) * Ls + LD / 2 - L2p / 2 - Ld1 / 2
    tot = sp.expand(sp.simplify(logR - logm + logw))
    check("(G6) log(Q^O_j/Q^G_j) - log m(theta_j) + log omega = 0 identically (symbolic s, r, theta, sigma, d, a)", tot == 0,
          f"residual {tot}")
    # tamper: omega with g -> g + 1/7 must leave a residual
    tt = sp.expand(sp.simplify(logR - logm + logw + Ls / 7))
    check("(G6) tamper g -> g + 1/7 in omega leaves a nonzero residual", tt != 0, f"{tt}")
    # tamper: Gauss with k = a instead of a + 1 (i.e. forgetting the missing factor) breaks G2
    check("(G2) tamper: u + (a-1)/a != beta1' - al_j", sp.simplify(u + (a - 1) / a - (b1 - al)) != 0)


def numeric():
    mp.mp.dps = 100
    worst = mp.mpf(0)
    for k in range(2, 10):
        for uu in (mp.mpf("0.3141592653"), mp.mpf("1.7320508"), mp.mpf("-0.45") + mp.mpf("0.2") * 1j):
            lhs = mp.fprod(mp.gamma(uu + mp.mpf(j) / k) for j in range(k))
            rhs = (2 * mp.pi) ** (mp.mpf(k - 1) / 2) * mp.power(k, mp.mpf(1) / 2 - k * uu) * mp.gamma(k * uu)
            worst = max(worst, abs(lhs / rhs - 1))
    check("(GM) Gauss multiplication at dps 100, k = 2..9, three u (one complex)", worst < mp.mpf(10) ** -95, mp.nstr(worst, 3))
    mp.mp.dps = 60
    MOM = {"A": (sp.Rational(1, 3), sp.Rational(3, 5)), "B": (sp.Rational(9, 10), sp.Rational(1, 5))}
    for sol, ai in (("Sol3", 3), ("Sol1", 3), ("Sol3", 1), ("Sol1", 1)):
        pt = SL.point(sol, ai)
        for nm, (l0, l1) in MOM.items():
            th, r, frozen = SL.exponents(pt, l0, l1)
            QO, _ = SL.meijer_Q0_ode(pt, th, r)
            QG, _ = SL.meijer_Q0_gamma(pt, th, frozen)
            w = mp.mpf(sp.N(SL.omega_sealed(pt, r), 70))
            dev = max(abs(QO[h] * w / (mp.mpf(sp.N(SL.m_sym(pt, r, th[h]), 70)) * QG[h]) - 1) for h in range(len(th)))
            check(f"lib tie {sol} a={ai} {nm}: omega Q^O(0)/(m Q^G(0)) = 1 (dps 60)", dev < mp.mpf(10) ** -50, mp.nstr(dev, 3))


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] != "GAUSS":
        raise SystemExit("usage: spec1_s2c_gauss.py GAUSS")
    SL.banner(__file__)
    symbolic()
    numeric()
    SL.imported("imported(end)")
    print("GAUSS DONE: " + ("ALL PASS" if not FAIL else f"FAILURES {FAIL}"), flush=True)
