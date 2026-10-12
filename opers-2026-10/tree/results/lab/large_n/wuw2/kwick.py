"""WUW2 engine: WUW1's free-field routines with contraction normalisation KAPPA (d Phi_c d Phi_c ~ KAPPA/(z-w)^2), plus sparse exact elimination."""
from fractions import Fraction as F
from math import factorial
from itertools import combinations
import wick as Wk
KAPPA = F(1)
def set_kappa(k):
    global KAPPA; KAPPA = F(k)
A, B = 1, 2
def weight(m): return sum(q + 1 for (_, q) in m)
def monomials(w, comps=(A, B)):
    letters = [(c, q) for c in comps for q in range(w)]
    out = []
    def rec(start, rem, acc):
        if rem == 0: out.append(tuple(sorted(acc))); return
        for i in range(start, len(letters)):
            c, q = letters[i]
            if q + 1 <= rem: rec(i, rem - q - 1, acc + [(c, q)])
    rec(0, w, []); return sorted(set(out))
def dpoly(G): return Wk.dfield(G)
def add(X, Y, s=1): return Wk.fadd(X, Y, s)
def scale(X, s): return Wk.fscale(X, s)
def mul(X, Y):
    out = {}
    for m1, c1 in X.items():
        for m2, c2 in Y.items():
            m = tuple(sorted(m1 + m2)); out[m] = out.get(m, 0) + c1 * c2
    return {m: c for m, c in out.items() if c}
def ope(X, Y, poles=None):
    """Wick OPE with KAPPA per contraction"""
    out = {}
    for mA, cA in X.items():
        for mB, cB in Y.items():
            for pairs, remA, remB in Wk._matchings(mA, mB):
                P = 0; fac = F(cA * cB)
                for (i, p), (_, q) in pairs:
                    fac *= KAPPA * (-1)**p * factorial(p + q + 1); P += p + q + 2
                for T in range(P + 1):
                    if poles is not None and (P - T) not in poles: continue
                    for sh, cf in Wk._taylor(remA, T):
                        mono = tuple(sorted(sh + list(remB))); r = P - T
                        d = out.setdefault(r, {}); d[mono] = d.get(mono, 0) + fac * cf
    return {r: {m: c for m, c in f.items() if c} for r, f in out.items()}
def residue_exp(G, alpha):
    out = {}
    for m, c in G.items():
        n = len(m)
        for k in range(1, n + 1):
            for S in combinations(range(n), k):
                fac = F(c); P = 0
                for i in S:
                    cc, q = m[i]; fac *= KAPPA * alpha[cc] * ((-1) ** q * factorial(q)); P += q + 1
                rem = [m[i] for i in range(n) if i not in S]
                for sh, cf in Wk._taylor(rem, P - 1):
                    nm = tuple(sorted(sh)); out[nm] = out.get(nm, 0) + fac * cf
    return {k: v for k, v in out.items() if v}
def D_exp(Qm, alpha):
    out = dict(dpoly({Qm: F(1)}))
    for cc in (A, B):
        if alpha[cc]: nm = tuple(sorted(Qm + ((cc, 0),))); out[nm] = out.get(nm, 0) + alpha[cc]
    return {k: v for k, v in out.items() if v}
def schur(alpha, comp, n):
    gen = {k: {((comp, k - 1),): F(alpha) / factorial(k)} for k in range(1, n + 1)}
    e = [{(): F(1)}]
    for m in range(1, n + 1):
        acc = {}
        for k in range(1, m + 1): acc = add(acc, scale(mul(gen[k], e[m - k]), F(k)))
        e.append(scale(acc, F(1, m)))
    return e[n]
def field_screen(G, alpha_c, comp):
    out = {}
    for m, c in G.items():
        n = len(m); idxs = [i for i in range(n) if m[i][0] == comp]
        for k in range(1, len(idxs) + 1):
            for S in combinations(idxs, k):
                fac = F(c); P = 0
                for i in S:
                    q = m[i][1]; fac *= -alpha_c * KAPPA * factorial(q); P += q + 1
                rem = tuple(m[i] for i in range(n) if i not in S)
                out = add(out, scale(mul(schur(alpha_c, comp, P - 1), {tuple(sorted(rem)): F(1)}), fac))
    return out
def current(a, comp):
    """Res e^(a Phi)(z) e^(-a Phi)(w) = S_(a^2 KAPPA - 1)(a d^k Phi/k!)"""
    n = F(a) ** 2 * KAPPA - 1; assert n.denominator == 1
    return schur(F(a), comp, int(n))
# ---------------- sparse exact echelon
class Ech:
    def __init__(s): s.rows = {}
    def reduce(s, v):
        """FULL reduction (every pivot key eliminated, largest first): a canonical, linear normal form modulo the span.
        (v1 eliminated only the leading key: correct for membership, NOT linear -- found by the kappa = 1 control; disclosed.)"""
        v = dict(v)
        while True:
            ks = [k for k in v if k in s.rows]
            if not ks: return v
            p = max(ks); c = v[p]; r = s.rows[p]
            for k, x in r.items():
                nv = v.get(k, 0) - c * x
                if nv: v[k] = nv
                else: v.pop(k, None)
    def insert(s, v):
        v = s.reduce(v)
        if not v: return False
        p = max(v); c = v[p]; s.rows[p] = {k: x / c for k, x in v.items()}; return True
def relations(vectors, ech=None):
    """basis of {x : sum x_i v_i in span(ech)} (dense over the reduced keys)"""
    red = [ech.reduce(v) if ech else dict(v) for v in vectors]
    keys = sorted({k for v in red for k in v}); idx = {k: i for i, k in enumerate(keys)}
    n = len(vectors)
    rows = [[F(0)] * n for _ in keys]
    for j, v in enumerate(red):
        for k, x in v.items(): rows[idx[k]][j] += x
    # rref
    M = rows; piv = []; r = 0
    for c in range(n):
        pr = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if pr is None: continue
        M[r], M[pr] = M[pr], M[r]; inv = 1 / M[r][c]; M[r] = [x * inv for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c]; M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        piv.append(c); r += 1
    free = [j for j in range(n) if j not in piv]; basis = []
    for f in free:
        x = [F(0)] * n; x[f] = F(1)
        for ri, pc in enumerate(piv): x[pc] = -M[ri][f]
        basis.append(x)
    return basis
def d_ech(W):
    e = Ech()
    for m in monomials(W - 1): e.insert(dpoly({m: F(1)}))
    return e
def commutant(W, alphas):
    """weight-W fields G with Res(G, alpha) = D_alpha(Q) for all alpha; returns list of fields (mod nothing)"""
    gm = monomials(W); qm = monomials(W - 2)
    vecs = []
    for g in gm:
        v = {}
        for ai, al in enumerate(alphas):
            for mono, x in residue_exp({g: F(1)}, al).items(): v[(ai, mono)] = v.get((ai, mono), 0) + x
        vecs.append(v)
    for ai, al in enumerate(alphas):
        for q in qm:
            vecs.append({(ai, mono): -x for mono, x in D_exp(q, al).items()})
    rel = relations(vecs)
    return [{gm[j]: x[j] for j in range(len(gm)) if x[j] != 0} for x in rel if any(x[j] != 0 for j in range(len(gm)))]
def reps_mod_d(fields, W):
    e = d_ech(W); out = []
    for f in fields:
        if e.insert(f): out.append(f)
    return out
