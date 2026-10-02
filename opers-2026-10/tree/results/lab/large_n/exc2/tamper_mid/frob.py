"""EXC2: local (Frobenius) analysis of a class-U operator with an apparent singularity at xi = z.

Operator (vartheta = xi d/dxi):
    L = P(vartheta) + xi Lm(vartheta) + Esign * E xi^b + xi sum_j R_j(xi) vartheta^j,   R_j = sum_m r[j][m]/(xi - z)^m.
Local variable d = xi - z.  Acting on d^e:   vartheta d^e = e (z d^(e-1) + d^e).
L d^e = sum_k l_k(e) d^(e+k),  k >= -order.   Frobenius: sum_k l_k(e + j - k) c_(j-k) = 0.
No-logarithm conditions for integer exponents e_1 < ... < e_d: start at e_1 with c_0 = 1, introduce a free constant at
every later exponent, and require the obstruction at each resonance to vanish identically in those constants and in E.
All exact (sympy)."""
import sys
import sympy as sp

sys.dont_write_bytecode = True
e, z, E, dlt = sp.symbols('e z E delta')


def theta_on(series):
    """series: dict k -> g_k(e) meaning sum g_k d^(e+k);  apply vartheta = (z + d) d/dd."""
    out = {}
    for k, g in series.items():
        f = sp.expand(g * (e + k))
        out[k - 1] = sp.expand(out.get(k - 1, 0) + z * f)
        out[k] = sp.expand(out.get(k, 0) + f)
    return {k: v for k, v in out.items() if v != 0}


def mul_series(series, coeffs):
    """multiply by a Laurent series in d: coeffs dict power -> coefficient"""
    out = {}
    for k, g in series.items():
        for p, c in coeffs.items():
            out[k + p] = sp.expand(out.get(k + p, 0) + g * c)
    return {k: v for k, v in out.items() if v != 0}


def add_series(*ss):
    out = {}
    for s in ss:
        for k, v in s.items():
            out[k] = sp.expand(out.get(k, 0) + v)
    return {k: v for k, v in out.items() if v != 0}


def trunc(series, K):
    return {k: v for k, v in series.items() if k <= K}


def local_operator(Pcoef, Lcoef, b, Esign, R, K):
    """Pcoef: list, P = sum Pcoef[j] vartheta^j; Lcoef likewise; R: dict j -> dict m -> coefficient of xi/(xi-z)^m vartheta^j.
    Returns dict k -> l_k(e), k <= K."""
    order = len(Pcoef) - 1
    powers = [{0: sp.Integer(1)}]
    for j in range(order):
        powers.append(trunc(theta_on(powers[-1]), K))
    xi = {0: z, 1: sp.Integer(1)}
    out = {}
    for j, c in enumerate(Pcoef):
        out = add_series(out, {k: c * v for k, v in powers[j].items()})
    for j, c in enumerate(Lcoef):
        out = add_series(out, mul_series({k: c * v for k, v in powers[j].items()}, xi))
    # E xi^b = E z^b (1 + d/z)^b ; write E0 = E z^b
    xib = {i: sp.binomial(b, i) / z**i for i in range(0, K + order + 1)}
    out = add_series(out, {k: Esign * E * v for k, v in xib.items()})
    for j, rj in R.items():
        for m, c in rj.items():
            out = add_series(out, mul_series(mul_series(powers[j], xi), {-m: c}))
    return trunc(out, K)


def nolog_conditions(lk, exps, order):
    """lk: dict k -> l_k(e) with k >= -order.  exps sorted integers.  Returns (indicial polynomial, list of conditions).
    Each condition is an expression that must vanish identically in E and in the free constants kap_i."""
    e1 = exps[0]
    span = exps[-1] - e1
    lead = lk.get(-order, 0)
    c = {0: sp.Integer(1)}
    kap = {}
    conds = []
    for j in range(1, span + 1):
        rhs = 0
        for kk in range(1, j + 1):
            lcoef = lk.get(kk - order, 0)
            if lcoef == 0 or (j - kk) not in c:
                continue
            rhs += lcoef.subs(e, e1 + j - kk) * c[j - kk]
        rhs = sp.expand(rhs)
        if (e1 + j) in exps:
            conds.append((e1 + j, rhs))
            kap[j] = sp.Symbol('kap%d' % j)
            c[j] = kap[j]
        else:
            den = lead.subs(e, e1 + j)
            c[j] = sp.cancel(-rhs / den)
    return lead, conds, list(kap.values())


def split_conditions(conds, kaps, extra_syms=()):
    """coefficients of each condition w.r.t. E and the free constants -> list of polynomial equations (numerators)."""
    eqs = []
    for pos, expr in conds:
        num = sp.numer(sp.together(expr))
        P = sp.Poly(sp.expand(num), E, *kaps)
        for mon, cf in P.terms():
            if cf != 0:
                eqs.append((pos, mon, sp.factor(cf)))
    return eqs
