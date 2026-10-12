"""VIR4-R -- own theta-form WKB engine (extends VIR3-R's BLZ-calibrated engine).  Equation prod_j (theta - lambda_j) psi = x^a (x^b - E) psi, u = log x,
v = e^(b u), E = -E': Q = v^A (v + E') with A = a/b.  psi = exp(int S du), S = sum_k Lambda^(1-k) s_k; the factors commute, so
F = sum_m c_m Y_(n-m), Y_0 = 1, Y_(r+1) = D Y_r + S Y_r, D = d/du = b v d/dv, c_m = coefficients of prod_j (t - lambda_j) = sum_m c_m t^(n-m).
s_0^n = U; s_k = -[F]_(Lambda^(n-k)) |_(s_k = 0) / (n s_0^(n-1)).  Expressions: {(alpha, beta): poly}, poly = {(i, j): Fraction} in (l0, l1) (plus a
generic parameter for the calibration).  Integrals: int du v^alpha (v + E')^beta = (1/|b|) E'^(alpha+beta) G(alpha) G(-alpha-beta)/G(-beta) (continuation).
Module."""
from fractions import Fraction as F
import sympy as sp
# ---- polynomials in (l0, l1) as dicts
def padd(p, q, s=1):
    out = dict(p)
    for k, v in q.items():
        out[k] = out.get(k, 0) + s * v
        if out[k] == 0: del out[k]
    return out
def pmul(p, q):
    out = {}
    for (i1, j1), a in p.items():
        for (i2, j2), b in q.items():
            k = (i1 + i2, j1 + j2); out[k] = out.get(k, 0) + a * b
    return {k: v for k, v in out.items() if v != 0}
def pscale(p, c): return {k: v * c for k, v in p.items() if v * c != 0}
# ---- expressions {(alpha, beta): poly}
def eadd(E1, E2, s=1):
    out = {k: dict(v) for k, v in E1.items()}
    for k, v in E2.items():
        out[k] = padd(out.get(k, {}), v, s)
        if not out[k]: del out[k]
    return out
def emul(E1, E2):
    out = {}
    for (a1, b1), p1 in E1.items():
        for (a2, b2), p2 in E2.items():
            k = (a1 + a2, b1 + b2); out[k] = padd(out.get(k, {}), pmul(p1, p2))
            if not out[k]: del out[k]
    return out
def escale(E1, c): return {k: pscale(v, c) for k, v in E1.items() if pscale(v, c)}
def eD(E1, bb):
    """d/du = b v d/dv on v^a (v+E')^c: b[a v^a (v+E')^c + c v^(a+1) (v+E')^(c-1)]"""
    out = {}
    for (a, c), p in E1.items():
        if a != 0: out = eadd(out, {(a, c): pscale(p, bb * a)})
        if c != 0: out = eadd(out, {(a + 1, c - 1): pscale(p, bb * c)})
    return out
def epoly_div_mono(E1, mono):   # divide by v^a0 (v+E')^c0
    a0, c0 = mono; return {(a - a0, c - c0): p for (a, c), p in E1.items()}
# ---- series in Lambda: {power: expr}
def smul(S1, S2, pmin):
    out = {}
    for p1, e1 in S1.items():
        for p2, e2 in S2.items():
            if p1 + p2 < pmin: continue
            out[p1 + p2] = eadd(out.get(p1 + p2, {}), emul(e1, e2))
    return {k: v for k, v in out.items() if v}
def sD(S1, bb): return {k: eD(v, bb) for k, v in S1.items() if eD(v, bb)}
def sadd(S1, S2, s=1):
    out = dict(S1)
    for k, v in S2.items():
        out[k] = eadd(out.get(k, {}), v, s)
        if not out[k]: del out[k]
    return out
def wkb(cm, n, A, bb, K, c0=F(1)):
    """cm: list of n+1 polys (c_0 = 1, ..., c_n); potential U = v^A (v + E'); returns [s_0..s_K] as exprs"""
    s0 = {(F(A, n), F(1, n)): {(0, 0): F(c0)}}          # c0^n = sigma (overall sign of the potential)
    s = [s0]
    den_mono = (F(A, n) * (n - 1), F(1, n) * (n - 1))
    for k in range(1, K + 1):
        pmin = n - k
        S = {1 - j: s[j] for j in range(k)}
        Y = [{0: {(F(0), F(0)): {(0, 0): F(1)}}}]
        for r in range(n):
            keep = pmin - (n - (r + 1))          # Y_(r+1) needs powers >= pmin - (n - r - 1): each later factor raises the power by <= 1
            Yn = sadd(sD(Y[-1], bb), smul(S, Y[-1], keep))
            Y.append({p: e for p, e in Yn.items() if p >= keep})
        Fk = {}
        for m in range(n + 1):
            if not cm[m]: continue
            e = Y[n - m].get(pmin, {})
            Fk = eadd(Fk, {kk: pmul(v, cm[m]) for kk, v in e.items()})
        sk = escale(epoly_div_mono(Fk, den_mono), F(-1, n) / F(c0) ** (n - 1))
        s.append(sk)
    return s
def integrate(expr, bb):
    """int du of expr = sum poly * v^a (v+E')^c -> (E' power, gamma prefactor at reference, poly) with all terms reduced to a common reference"""
    if not expr: return None
    items = sorted(expr.items())
    tot_pow = {a + c for (a, c) in expr}; assert len(tot_pow) == 1, tot_pow
    (a0, c0), _ = items[0]
    def poch(x, m):
        r = F(1)
        if m >= 0:
            for i in range(m): r *= (x + i)
        else:
            for i in range(1, -m + 1): r /= (x - i)
        return r
    out = {}
    for (a, c), p in items:
        m = a - a0; assert m.denominator == 1; m = int(m)
        # G(a)/G(a0) = (a0)_m ; G(-c0)/G(-c) with -c = -c0 + m -> 1/(-c0)_m ; G(-a-c) common
        ratio = poch(a0, m) / poch(-c0, m)
        out = padd(out, pscale(p, ratio))
    pref = sp.Rational(1, abs(bb)) * sp.gamma(sp.Rational(a0.numerator, a0.denominator)) * sp.gamma(sp.Rational(-(a0 + c0).numerator, (a0 + c0).denominator)) / sp.gamma(sp.Rational((-c0).numerator, (-c0).denominator))
    return (a0 + c0), pref, out
def charpoly_coeffs(lams, n):
    """lams: list of polys (each lambda_j as poly in l0, l1); returns c_0..c_n of prod (t - lambda_j)"""
    c = [{(0, 0): F(1)}]
    for lam in lams:
        new = [dict(x) for x in c] + [{}]
        for i in range(len(c)):
            new[i + 1] = padd(new[i + 1], pmul(c[i], pscale(lam, -1)))
        c = new
    return c

def wkb2(c1, n1, c2, n2, A, bb, K, c0=F(1)):
    """two-sided theta-form: sum_m c1_m Y_(n1-m) = Q sum_m c2_m Y_(n2-m), Q = Lambda^(n1-n2) W, W = v^A (v + E');
    s_0^(n1-n2) = W (s_0 = c0 W^(1/(n1-n2)), c0^(n1-n2) = sign); s_k = -[F1 - Q F2]_(Lambda^(n1-k)) / ((n1 - n2) s_0^(n1-1)).
    With n2 = 0, c2 = [1] this is wkb()."""
    D = n1 - n2
    s0 = {(F(A, D), F(1, D)): {(0, 0): F(c0)}}
    s = [s0]
    den_mono = (F(A, D) * (n1 - 1), F(1, D) * (n1 - 1))
    W = {(F(A), F(1)): {(0, 0): F(1)}}
    nmax = max(n1, n2)
    for k in range(1, K + 1):
        S = {1 - j: s[j] for j in range(k)}
        Y = [{0: {(F(0), F(0)): {(0, 0): F(1)}}}]
        for r in range(nmax):
            keep = (r + 1) - k
            Yn = sadd(sD(Y[-1], bb), smul(S, Y[-1], keep))
            Y.append({p: e for p, e in Yn.items() if p >= keep})
        Fk = {}
        for m in range(n1 + 1):
            if c1[m]:
                e = Y[n1 - m].get(n1 - k, {})
                Fk = eadd(Fk, {kk: pmul(v, c1[m]) for kk, v in e.items()})
        for m in range(n2 + 1):
            if c2[m]:
                e = Y[n2 - m].get(n2 - k, {})
                Fk = eadd(Fk, emul(W, {kk: pmul(v, c2[m]) for kk, v in e.items()}), -1)
        sk = escale(epoly_div_mono(Fk, den_mono), F(-1, D) / F(c0) ** (n1 - 1))
        s.append(sk)
    return s

def integrate_reg(expr, bb):
    """the M-deformation limit of integrate(): when every term has alpha a non-positive integer (simple pole of Gamma(alpha)), the uniform 1/eps cancels
    in any monic normalisation and the limit is given by replacing Gamma(alpha = -m) with its residue (-1)^m/m!.  Returns (power, poly) with the common
    finite factor Gamma(-alpha-beta) divided out (constant within one s_k); terms without a pole are subleading and dropped (reported)."""
    import sympy as _sp
    from math import factorial as _fac
    if not expr: return None
    tot = {a + c for (a, c) in expr}; assert len(tot) == 1
    out = {}; nonpole = 0
    for (a, c), p in expr.items():
        if a.denominator == 1 and a <= 0:
            m = int(-a); res = _sp.Rational((-1) ** m, _fac(m))
            g = res / _sp.gamma(_sp.Rational((-c).numerator, (-c).denominator)) if not ((-c).denominator == 1 and -c <= 0) else 0
            if g == 0: continue
            # Gamma(-c) values for half-integers share sqrt(pi): factor it out via ratio to Gamma(1/2)
            gr = _sp.simplify(g * _sp.sqrt(_sp.pi)) if (-c).denominator == 2 else _sp.simplify(g)
            for kk, v in p.items(): out[kk] = out.get(kk, 0) + _sp.Rational(v.numerator, v.denominator) * gr
        else:
            nonpole += 1
    return tot.pop(), {k: v for k, v in out.items() if v != 0}, nonpole
