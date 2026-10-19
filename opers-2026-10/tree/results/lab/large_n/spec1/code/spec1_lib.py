"""SPEC1 shared library.  NEW FILE.
Exact parameters of the S0 family (deg L = 1, positive integer a), the multiplier m, the sealed constant omega (SPEC1_PREREG sec. 4),
the beta formulas, and the E = 0 Meijer-G connection coefficients of both forms (S2c, derived in the PREREG's plan; see spec1_s2.py).
nsimplify / lambdify disabled."""
import hashlib
import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, os.pardir))
sys.path.insert(0, os.path.join(ROOT, "deps_snapshot"))

import sympy as sp  # noqa: E402
import mpmath as mp  # noqa: E402


def _no(*a, **k):
    raise RuntimeError("nsimplify / lambdify are forbidden in SPEC1")


sp.nsimplify = _no
sp.lambdify = _no


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def banner(script):
    print(f"script {os.path.abspath(script)} sha256 {sha(script)}", flush=True)
    print(f"python {sys.version.split()[0]} sympy {sp.__version__} mpmath {mp.__version__}", flush=True)
    imported()


def imported(tag="imported"):
    seen = set()
    for _, mod in sorted(sys.modules.items()):
        f = getattr(mod, "__file__", None)
        if f and os.path.abspath(f).startswith(ROOT) and f.endswith(".py") and os.path.abspath(f) not in seen:
            seen.add(os.path.abspath(f))
            print(f"{tag} {os.path.relpath(os.path.abspath(f), ROOT)} sha256 {sha(f)}", flush=True)


def point(sol, a):
    """Exact parameters (Fractions) of the S0 point (sol in {'Sol3', 'Sol1'}, positive integer a)."""
    a = F(a)
    if sol == "Sol3":
        t = (3 * a + 2) / (a + 2)
        n, d = a + 4, 4
    else:
        t = 2 * a + 1
        n, d = (t + 3) / 2, 2
    sigma = n / a
    s = (n - 1) / 2
    c = sigma * (a + 1)
    return dict(sol=sol, a=a, t=t, n=n, d=d, D=d + int(a), sigma=sigma, s=s, c=c)


def exponents(pt, l0, l1):
    """(dynamical theta_j, r, frozen exponents) as sympy numbers (exact if l0, l1 are rational)."""
    s, sigma, a = sp.Rational(pt["s"]), sp.Rational(pt["sigma"]), int(pt["a"])
    if pt["sol"] == "Sol3":
        r = sp.Integer(0)
        th = [s - l0, s + l0, s - l1, s + l1]
    else:
        r = l0
        th = [s - l1, s + l1]
    base = s - r + sigma * (1 - a) / 2
    frozen = [base + sigma * j for j in range(a)]
    return th, r, frozen


def g_of(pt, r):
    return (sp.Rational(pt["s"]) - r) / sp.Rational(pt["sigma"]) + sp.Rational(1 - int(pt["a"]), 2)


def m_sym(pt, r, nu):
    sigma = sp.Rational(pt["sigma"])
    return sigma ** (nu / sigma) / sp.gamma(g_of(pt, r) - nu / sigma)


def omega_sealed(pt, r):
    """SPEC1_PREREG sec. 4 (SEALED): omega = sigma^(g - 1/2) sqrt((d + a) / (2 pi (d - 1)))."""
    sigma = sp.Rational(pt["sigma"])
    d, a = pt["d"], int(pt["a"])
    return sigma ** (g_of(pt, r) - sp.Rational(1, 2)) * sp.sqrt(sp.Rational(d + a, 2 * (d - 1)) / sp.pi)


def beta_gamma(pt, r):
    s, D, a = sp.Rational(pt["s"]), pt["D"], int(pt["a"])
    qG = sp.Rational(pt["c"]) / D
    return s - qG * (D - 1) / 2 - a * r / D


def beta_ode(pt, r):
    s, c, d = sp.Rational(pt["s"]), sp.Rational(pt["c"]), pt["d"]
    return s - c / 2 + r / (d - 1)


# ------------------------------------------------------------------ E = 0 closed forms (mpmath), S2c
def meijer_Q0_ode(pt, th, r):
    """Q^O_h(0) for every dynamical h.  At E = 0: P(s-mu) a_mu = (s - mu + c/2 - r) a_(mu-c), so the subdominant solution is
    G^{d,0}_{1,d}(c^(1-d) x^c | beta1; alpha_1..alpha_d), alpha = theta/c, beta1 = (s + c/2 - r)/c.  Frobenius coefficient
    prod_{i != h} Gamma(alpha_i - alpha_h)/Gamma(beta1 - alpha_h) * c^((1-d) alpha_h); asymptotics (nu = d - 1)
    (2 pi)^((nu-1)/2) nu^(-1/2) exp(-nu w^(1/nu)) w^vt, vt = ((1 - nu)/2 + sum alpha - beta1)/nu, w = c^(1-d) x^c."""
    c = mp.mpf(sp.Rational(pt["c"]).p) / sp.Rational(pt["c"]).q
    d = pt["d"]
    s = mp.mpf(sp.Rational(pt["s"]).p) / sp.Rational(pt["s"]).q
    al = [mp.mpf(sp.N(x, mp.mp.dps + 10)) / c for x in th]
    rr = mp.mpf(sp.N(r, mp.mp.dps + 10))
    b1 = (s + c / 2 - rr) / c
    nu = d - 1
    vt = ((1 - mp.mpf(nu)) / 2 + sum(al) - b1) / nu
    norm = (2 * mp.pi) ** (mp.mpf(nu - 1) / 2) / mp.sqrt(nu) * c ** ((1 - d) * vt)
    out = []
    for h in range(len(al)):
        v = mp.mpf(1)
        for i in range(len(al)):
            if i != h:
                v *= mp.gamma(al[i] - al[h])
        out.append(v / mp.gamma(b1 - al[h]) * c ** ((1 - d) * al[h]) / norm)
    return out, c * vt


def meijer_Q0_gamma(pt, th, frozen):
    """Q^Gamma_h(0) for the dynamical h.  At E = 0: Pi(s - mu) c_mu = c_(mu - c) with Pi of degree D, so the subdominant solution is
    G^{D,0}_{0,D}(c^(-D) x^c | ; eps_1..eps_D), eps = (dynamical + frozen exponents)/c; asymptotics nu = D,
    vt = ((1 - D)/2 + sum eps)/D."""
    c = mp.mpf(sp.Rational(pt["c"]).p) / sp.Rational(pt["c"]).q
    D = pt["D"]
    ex = [mp.mpf(sp.N(x, mp.mp.dps + 10)) for x in list(th) + list(frozen)]
    ep = [x / c for x in ex]
    vt = ((1 - mp.mpf(D)) / 2 + sum(ep)) / D
    norm = (2 * mp.pi) ** (mp.mpf(D - 1) / 2) / mp.sqrt(D) * c ** (-D * vt)
    out = []
    for h in range(len(th)):
        v = mp.mpf(1)
        for i in range(len(ep)):
            if i != h:
                v *= mp.gamma(ep[i] - ep[h])
        out.append(v * c ** (-D * ep[h]) / norm)
    return out, c * vt
