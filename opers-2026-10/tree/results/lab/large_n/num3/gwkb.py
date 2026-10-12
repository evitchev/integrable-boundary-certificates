"""SOL12 -- own general theta-symbol WKB.  Equation G(s - theta) psi = x^n (x^(nM) - E) psi, u = log x, v = x^(nM), E = -E', W = v^A (v + E'), A = 1/M,
d/du = bb v d/dv with bb = nM.  G(sigma) = sigma^D exp(sum_r gamma_r sigma^-r) asymptotically (sigma = s + Sigma), built from blocks
R_a((sigma - shift)/delta) delta^a, R_a(u) = Gamma(u + (a+1)/2)/Gamma(u - (a-1)/2) (Stirling/Bernoulli).
psi = exp(int S du), S = -Sigma; f(theta + S)1 = sum_m c_m f^(m)(S), c_m = [a^m] exp(sum_(j>=2) a^j S^(j-1)/j!), f^(m)(S) = (-1)^m G^(m)(s + Sigma).
Log form: D log(1+T) + Gt(Sigma) + log(1 + sum_(m>=2) (-1)^m rho_m c_m) = 0, Sigma = Lambda s0 (1 + T), s0^D = W, rho_m = G^(m)/G.
Expressions/polys as in thwkb.  Returns s_k = -s0 t_k (s_0 = -s0).  Module."""
from fractions import Fraction as F
import sympy as sp
from thwkb import padd, pmul, pscale, eadd, emul, escale, eD, sadd, smul, sD

def ppow_list(w, R):
    out = [{(0, 0): F(1)}]
    for _ in range(R): out.append(pmul(out[-1], w))
    return out
def binom(x, p):
    r = F(1)
    for i in range(p): r = r * (x - i) / (i + 1)
    return r
def bern(nn, x):
    v = sp.bernoulli(nn, sp.Rational(x.numerator, x.denominator)); return F(int(sp.fraction(v)[0]), int(sp.fraction(v)[1]))
def gammas(blocks, s, R):
    """blocks: list of (a, shift poly, delta).  Returns D, [gamma_1..gamma_R] (polys) with log G(s + Sigma) = D log Sigma + sum gamma_r Sigma^-r"""
    D = sum(F(b[0]) for b in blocks); g = [dict() for _ in range(R + 1)]
    for a, shift, delta in blocks:
        a = F(a); delta = F(delta)
        w = padd({(0, 0): F(s)}, shift, -1); wp = ppow_list(w, R)
        A_, B_ = (a + 1) / 2, -(a - 1) / 2
        bj = [None] + [F((-1) ** (j + 1)) * (bern(j + 1, A_) - bern(j + 1, B_)) / (j * (j + 1)) for j in range(1, R + 1)]
        for r in range(1, R + 1):
            g[r] = padd(g[r], pscale(wp[r], a * F((-1) ** (r + 1), r)))
            for j in range(1, r + 1):
                if bj[j] == 0: continue
                g[r] = padd(g[r], pscale(wp[r - j], bj[j] * delta ** j * binom(-j, r - j)))
    # the blocks come in +- pairs, so G is even in l0 and l1: work in (l0^2, l1^2)
    # (when G is not even, e.g. the BLZ calibration with shifts +-(l + 1/2), the original variables are kept)
    even = all(i % 2 == 0 and j % 2 == 0 for r in range(1, R + 1) for (i, j) in g[r])
    if even:
        for r in range(1, R + 1): g[r] = {(i // 2, j // 2): v for (i, j), v in g[r].items()}
    return D, g, even
def rhos(D, g, R):
    """rho_m(z) = G^(m)/G as z-series (z = 1/Sigma), m = 0..R, powers <= R: dict r -> poly"""
    lam = {1: {(0, 0): F(D)}}
    for r in range(1, R):
        if g[r]: lam[r + 1] = pscale(g[r], -r)
    rho = [{0: {(0, 0): F(1)}}]
    for m in range(R):
        cur = rho[-1]; nxt = {}
        for r, p in cur.items():
            if r + 1 <= R and r != 0: nxt[r + 1] = padd(nxt.get(r + 1, {}), pscale(p, -r))
            for rr, q in lam.items():
                if r + rr <= R: nxt[r + rr] = padd(nxt.get(r + rr, {}), pmul(p, q))
        rho.append({k: v for k, v in nxt.items() if v})
    return rho
def strunc(S1, pmin): return {k: v for k, v in S1.items() if k >= pmin and v}
def sscale_expr(S1, e): return {k: emul(v, e) for k, v in S1.items() if emul(v, e)}
def spolymul(S1, p): return {k: {kk: pmul(vv, p) for kk, vv in v.items() if pmul(vv, p)} for k, v in S1.items()}
def slog1p(X, pmin, K):
    out = {}; P = {0: {(F(0), F(0)): {(0, 0): F(1)}}}
    for q in range(1, K + 1):
        P = strunc(smul(P, X, pmin), pmin)
        if not P: break
        out = sadd(out, {k: escale(v, F((-1) ** (q + 1), q)) for k, v in P.items()})
    return out
def wkb(blocks, s, M, K, verbose=False):
    n = sum(F(b[0]) for b in blocks); M = F(M); A = 1 / M; bb = n * M
    R = 2 * K + 2
    D, g, even = gammas(blocks, s, R); rho = rhos(D, g, R)
    s0 = {(A / D, F(1) / D): {(0, 0): F(1)}}
    one = {(F(0), F(0)): {(0, 0): F(1)}}
    t = [None]
    for k in range(1, K + 1):
        pmin = -k
        T = {-j: t[j] for j in range(1, k)}
        Tp = [{0: one}]
        for p in range(1, k + 1): Tp.append(strunc(smul(Tp[-1], T, pmin), pmin))
        tot = {}
        # D log(1+T)
        tot = sadd(tot, {kk: escale(v, D) for kk, v in slog1p(T, pmin, k).items()})
        # (1+T)^(-r) and Sigma^(-r) = Lambda^-r s0^-r (1+T)^-r
        def one_plus_T_pow(r):
            out = {}
            for p in range(0, k + 1):
                c = binom(F(-r), p)
                if c and Tp[p]: out = sadd(out, {kk: escale(v, c) for kk, v in Tp[p].items()})
            return out
        Sinv = {}
        for r in range(1, R + 1):
            if r > k + R: break
            base = one_plus_T_pow(r); mono = {(-A / D * r, -F(r) / D): {(0, 0): F(1)}}
            Sinv[r] = strunc({kk - r: emul(v, mono) for kk, v in base.items()}, pmin - (R))
        # Gt(Sigma) = sum gamma_r Sigma^-r
        for r in range(1, k + 1):
            if g[r]: tot = sadd(tot, strunc(spolymul(Sinv[r], g[r]), pmin))
        # c_m from S = -Sigma
        Sig = {1: s0}
        for j in range(1, k): Sig[1 - j] = emul(s0, t[j])
        e = {}; Dj = {kk: escale(v, F(-1)) for kk, v in Sig.items()}
        fact = 1
        for j in range(2, 2 * k + 1):
            Dj = sD(Dj, bb); fact *= j
            e[j] = {kk: escale(v, F(1, fact)) for kk, v in Dj.items() if kk >= pmin - j}   # S^(j-1)/j!
        c = [{0: one}, {}]
        X = {}
        for m in range(2, 2 * k + 1):
            cm = {}
            for j in range(2, m + 1):
                if c[m - j]: cm = sadd(cm, {kk: escale(v, F(j, m)) for kk, v in smul(e[j], c[m - j], pmin + m).items()})
            cm = strunc(cm, pmin + m)
            c.append(cm)
            if not cm: continue
            for r, q in rho[m].items():
                if r < m or r > k + m // 2 + 1: continue
                term = smul(spolymul(Sinv[r], q), cm, pmin)
                X = sadd(X, {kk: escale(v, F((-1) ** m)) for kk, v in term.items()})
        X = strunc(X, pmin)
        tot = sadd(tot, slog1p(X, pmin, k))
        tk = escale(tot.get(pmin, {}), F(-1) / D)
        assert all(kk >= pmin for kk in tot), 'positive powers?'
        bad = [kk for kk in tot if kk > pmin and kk <= 0 and tot[kk]]
        assert not bad, f'lower orders not cancelled at step {k}: {bad}'
        t.append(tk)
        if verbose: print(f'   step {k}: {len(tk)} terms', flush=True)
    def back(e): return e if not even else {kk: {(2 * i, 2 * j): v for (i, j), v in p.items()} for kk, p in e.items()}
    return [escale(s0, F(-1))] + [back(escale(emul(s0, t[k]), F(-1))) for k in range(1, K + 1)]
