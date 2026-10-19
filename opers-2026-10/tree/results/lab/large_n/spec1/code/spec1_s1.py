"""SPEC1 S1 (SPEC1_PREREG sec. 3): the general kernel.  NEW FILE.  Usage: spec1_s1.py S1
(a) m(mu + sigma)/m(mu) = s - mu - r - c/2 (c = sigma (a + 1)), symbolic in a, sigma, s, r, mu.  Tamper g + 1/7.
(b) m(mu)/m(mu - n) = G_a(s - mu - r) (n = sigma a) for general real a (argument identity of the Gamma ratio), and for integer
    a = 1..8 as the expanded frozen product prod_j (T - sigma (a-1)/2 + sigma j).  Tamper g + 1/7.
(c) frozen exponents base + sigma j (j < a) are roots of G_a(s - f - r) and poles of 1/Gamma(g - f/sigma) (argument -j), a = 1..8;
    the frozen Frobenius lattice (steps sigma a, sigma (a+1)) stays in base + sigma Z>=0 and misses the other frozen exponents.
    Tamper base + sigma/2.
(d) intertwining on monomials, L_O K x^mu = x^n K x^-n L_Gamma x^mu, coefficientwise, symbolic in mu, l0, l1, E, at every S3 point.
(e) quadrature (mpmath, 50 digits, log-tau parametrised Hankel contour): K x^mu = m(mu) x^mu at complex mu, and Lemma 2(i)
    K_(g-1) psi = K_g[(g - 1 - theta/sigma) psi] for psi = x^0.4 exp(-x^2), at Sol3 a = 1, 3 and Sol1 a = 3 (point A)."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import spec1_lib as SL  # noqa: E402
from spec1_lib import sp, mp, F  # noqa: E402

FAIL = []
POINTS = [("Sol3", 1), ("Sol3", 2), ("Sol3", 3), ("Sol1", 1), ("Sol1", 2), ("Sol1", 3)]


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""), flush=True)
    if not ok:
        FAIL.append(name)


def mratio(g, sig, nu1, nu2):
    """run 2: m(nu1)/m(nu2) = sigma^((nu1-nu2)/sigma) Gamma(g - nu2/sigma)/Gamma(g - nu1/sigma), the Gamma ratio reduced to a finite
    product (its argument difference must be an integer; sympy's gammasimp did not reduce it in run 1)."""
    A, B = g - nu2 / sig, g - nu1 / sig
    k = sp.simplify(A - B)
    assert k.is_integer, k
    k = int(k)
    gr = sp.prod([B + i for i in range(k)]) if k >= 0 else 1 / sp.prod([A + i for i in range(-k)])
    return sp.powsimp(sig ** ((nu1 - nu2) / sig), force=True) * gr


def part_a_b():
    a, sig, s, r, mu = sp.symbols("a sigma s r mu", positive=True)
    g = (s - r) / sig + (1 - a) / 2
    c = sig * (a + 1)
    m = lambda nu, gg=g: sig ** (nu / sig) / sp.gamma(gg - nu / sig)  # noqa: E731
    lhs = sp.simplify(mratio(g, sig, mu + sig, mu))
    check("(a) m(mu+sigma)/m(mu) = s - mu - r - c/2 (symbolic)", sp.simplify(lhs - (s - mu - r - c / 2)) == 0, f"{lhs}")
    gt = g + sp.Rational(1, 7)
    lhs_t = sp.simplify(mratio(gt, sig, mu + sig, mu))
    check("(a) tamper g + 1/7 fails", sp.simplify(lhs_t - (s - mu - r - c / 2)) != 0)
    # (b) general a: m(mu)/m(mu - n) = sigma^a Gamma(g - mu/sigma + a)/Gamma(g - mu/sigma); G_a(s-mu-r) = sigma^a Gamma(u+(a+1)/2)/Gamma(u-(a-1)/2)
    u = (s - mu - r) / sig
    n = sig * a
    A1 = sp.simplify((g - (mu - n) / sig) - (u + (a + 1) / 2))
    A2 = sp.simplify((g - mu / sig) - (u - (a - 1) / 2))
    pw = sp.simplify(sp.powsimp(sig ** (mu / sig) / sig ** ((mu - n) / sig), force=True) - sig ** a)
    check("(b) general a: Gamma arguments and the power of sigma match exactly", A1 == 0 and A2 == 0 and pw == 0,
          f"{A1}, {A2}, {pw}")
    A1t = sp.simplify((gt - (mu - n) / sig) - (u + (a + 1) / 2))
    check("(b) tamper g + 1/7 fails", A1t != 0)
    T = sp.Symbol("T")
    okb = True
    for ai in range(1, 9):
        Ga_poch = sig ** ai * sp.rf(T / sig - sp.Rational(ai - 1, 2), ai)
        Ga_prod = sp.prod([T - sig * sp.Rational(ai - 1, 2) + sig * j for j in range(ai)])
        okb &= sp.expand(Ga_poch - Ga_prod) == 0
        ratio = sp.expand(mratio(g.subs(a, ai), sig, mu, mu - sig * ai))
        okb &= sp.expand(ratio - Ga_prod.subs(T, s - mu - r)) == 0
    check("(b) integer a = 1..8: m(mu)/m(mu-n) = prod_j (T - sigma(a-1)/2 + sigma j) at T = s - mu - r (expanded)", okb)


def part_c():
    sig, s, r = sp.symbols("sigma s r", positive=True)
    T = sp.Symbol("T")
    ok = okt = True
    lat_ok = True
    for ai in range(1, 9):
        g = (s - r) / sig + sp.Rational(1 - ai, 2)
        base = s - r + sig * sp.Rational(1 - ai, 2)
        Ga = sp.prod([T - sig * sp.Rational(ai - 1, 2) + sig * j for j in range(ai)])
        for j in range(ai):
            f = base + sig * j
            ok &= sp.expand(Ga.subs(T, s - f - r)) == 0 and sp.simplify(g - f / sig) == -j
            ft = f + sig / 2
            okt &= sp.simplify(g - ft / sig) != -j and sp.expand(Ga.subs(T, s - ft - r)) != 0
        # lattice: base + sigma j + sigma (ai k1 + (ai+1) k2), (k1, k2) != (0, 0), never equals another frozen exponent
        steps = sorted({ai * k1 + (ai + 1) * k2 for k1 in range(0, 4) for k2 in range(0, 4)} - {0})
        lat_ok &= all(j + st >= ai for j in range(ai) for st in steps)
    check("(c) frozen exponents base + sigma j are roots of G_a(s - f - r) and poles of 1/Gamma (argument -j), a = 1..8", ok)
    check("(c) tamper base + sigma/2: neither a root nor a pole", okt)
    check("(c) frozen Frobenius lattice (steps sigma a, sigma(a+1)) misses the other frozen exponents, a = 1..8", lat_ok)


def part_d():
    mu, l0, l1, E = sp.symbols("mu l0 l1 E")
    ok = True
    for sol, ai in POINTS:
        pt = SL.point(sol, ai)
        s, sig, c, n = (sp.Rational(pt[k]) for k in ("s", "sigma", "c", "n"))
        T = sp.Symbol("T")
        if sol == "Sol3":
            P = (T ** 2 - l0 ** 2) * (T ** 2 - l1 ** 2)
            r = sp.Integer(0)
        else:
            P = T ** 2 - l1 ** 2
            r = l0
        Ga = sp.prod([T - r - sig * sp.Rational(ai - 1, 2) + sig * j for j in range(ai)])
        gg = SL.g_of(pt, r)
        # L_O K x^mu = m(mu)[P(s-mu) x^mu - (s - mu - r - c/2) x^(mu+c) - E x^(mu+n)]
        # x^n K x^-n L_G x^mu = m(mu-n) P(s-mu) Ga(s-mu) x^mu - m(mu+c-n) x^(mu+c) - E m(mu) x^(mu+n)
        e1 = mratio(gg, sig, mu - n, mu) * P.subs(T, s - mu) * Ga.subs(T, s - mu) - P.subs(T, s - mu)
        e2 = mratio(gg, sig, mu + c - n, mu) - (s - mu - r - c / 2)
        ok_pt = sp.simplify(sp.expand(e1)) == 0 and sp.simplify(sp.expand(e2)) == 0
        print(f"   (d) {sol} a={ai}: x^mu coefficient diff {sp.simplify(e1)}, x^(mu+c) diff {sp.simplify(e2)}", flush=True)
        ok &= ok_pt
    check("(d) L_O K = x^n K x^-n L_Gamma on monomials at all six S0 points (symbolic in mu, l0, l1, E)", ok)


def hankel(g, sigma, x, psi, U=mp.mpf(240)):
    """(1/2 pi i) int_H e^tau tau^(-g) psi(log X) dtau, with log X = log x + (log sigma + L)/sigma, L = log tau (explicit)."""
    lx = mp.log(x)
    ls = mp.log(sigma)

    def integrand(L):
        tau = mp.exp(L)
        return mp.exp(tau) * mp.exp((1 - g) * L) * psi(lx + (ls + L) / sigma)
    lo = mp.quad(lambda u: integrand(mp.log(u) - 1j * mp.pi) / u, [U, 50, 10, 1])        # u: U -> 1 (lower bank)
    ci = mp.quad(lambda ph: integrand(1j * ph) * 1j, [-mp.pi, 0, mp.pi])                    # |tau| = 1
    up = mp.quad(lambda u: integrand(mp.log(u) + 1j * mp.pi) / u, [1, 10, 50, U])         # 1 -> U (upper bank)
    return (lo + ci + up) / (2j * mp.pi)


def part_e():
    mp.mp.dps = 50
    for sol, ai in (("Sol3", 1), ("Sol3", 3), ("Sol1", 3)):
        pt = SL.point(sol, ai)
        th, r, _ = SL.exponents(pt, sp.Rational(1, 3), sp.Rational(3, 5))
        g = mp.mpf(sp.N(SL.g_of(pt, r), 60))
        sig = mp.mpf(sp.N(sp.Rational(pt["sigma"]), 60))
        x = mp.mpf("0.7")
        nu = mp.mpc("0.37", "0.2")
        mono = hankel(g, sig, x, lambda lX: mp.exp(nu * lX))
        pred = sig ** (nu / sig) / mp.gamma(g - nu / sig) * x ** nu
        e_m = abs(mono / pred - 1)
        # run 2: psi = x^0.4 exp(-x) (growth order 1 < sigma at every point; run 1's exp(-x^2) has order 2 > sigma = 5/3 at Sol1 a=3)
        psi = lambda lX: mp.exp(mp.mpf("0.4") * lX - mp.exp(lX))  # noqa: E731
        tpsi = lambda lX: (mp.mpf("0.4") - mp.exp(lX)) * psi(lX)  # theta psi  # noqa: E731
        lhs = hankel(g - 1, sig, x, psi)
        rhs = hankel(g, sig, x, lambda lX: (g - 1) * psi(lX) - tpsi(lX) / sig)
        e_i = abs(lhs / rhs - 1)
        check(f"(e) {sol} a={ai} (g = {SL.g_of(pt, r)}, sigma = {pt['sigma']}): K x^mu = m(mu) x^mu at mu = 0.37+0.2i",
              e_m < mp.mpf(10) ** -40, f"rel err {mp.nstr(e_m, 3)}")
        check(f"(e) {sol} a={ai}: Lemma 2(i) on psi = x^0.4 exp(-x)", e_i < mp.mpf(10) ** -40,
              f"rel err {mp.nstr(e_i, 3)}; value {mp.nstr(lhs, 12)}")
        lhs_t = hankel(g - 1, sig, x, lambda lX: psi(lX) * mp.exp(lX / 7))
        check(f"(e) {sol} a={ai}: tamper (psi -> psi x^(1/7) on one side) differs", abs(lhs_t / rhs - 1) > mp.mpf(10) ** -5)


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] != "S1":
        raise SystemExit("usage: spec1_s1.py S1")
    SL.banner(__file__)
    part_a_b()
    part_c()
    part_d()
    part_e()
    SL.imported("imported(end)")
    print("S1 DONE: " + ("ALL PASS" if not FAIL else f"FAILURES {FAIL}"), flush=True)
