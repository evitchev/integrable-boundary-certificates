"""VIR2-R -- object: the record's cylindrical charges I_3, I_5, I_7 (mixed_engine representatives at curve_point(t)), rewritten in two bosons (X, phi)
via VIR1's Heis_X (x) Vir_N decomposition and T_Y -> T_phi = 1/2 (dphi)^2 + rho d^2 phi (rho = p/2 - 1/p, c = N); tested against Fable's VIR2 screening
pairs V_+- = exp(+-a X + beta phi), and the commutant of {exp(p phi), V_+, V_-} at weights 4, 6, 8.  Sealed: SEAL_VIR2R.md efa0c7f5....
Own code only (onwick.py formal-N pair Wick engine; wick.py explicit Wick engine; exact Q(sqrt d) linear algebra).  Usage: vir2r.py [--nocomm]"""
import os, sys, json, hashlib, time
from fractions import Fraction as F
from math import factorial, isqrt
from itertools import combinations
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
NOCOMM = '--nocomm' in sys.argv
LOG = open(os.path.join(HERE, 'run_VIR2R.log'), 'w'); T0 = time.time()
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); LOG.write(s + '\n'); LOG.flush()
FAIL = []; RES = []
def control(name, ok):
    say(f'CONTROL {name}: {ok}')
    if not ok: FAIL.append(name)
def result(name, ok):
    say(f'SEALED TEST {name}: {ok}')
    if not ok: RES.append(name)
assert hashlib.sha256(open(os.path.join(HERE, 'SEAL_VIR2R.md'), 'rb').read()).hexdigest().startswith('efa0c7f5'), 'seal'
assert open(os.path.join(HERE, 'run_ctrl_e.log')).read().strip().endswith('exit 0'), 'control (e) must pass first'
import onwick as O, wick as Wk
XC, FC = -1, -2
# ======================= exact field Q(sqrt d) =======================
class QS:
    __slots__ = ('u', 'v'); d = F(0)
    def __init__(s, u, v=0): s.u = F(u); s.v = F(v)
    def __add__(s, o): o = o if isinstance(o, QS) else QS(o); return QS(s.u + o.u, s.v + o.v)
    __radd__ = __add__
    def __sub__(s, o): o = o if isinstance(o, QS) else QS(o); return QS(s.u - o.u, s.v - o.v)
    def __rsub__(s, o): return QS(o) - s
    def __neg__(s): return QS(-s.u, -s.v)
    def __mul__(s, o): o = o if isinstance(o, QS) else QS(o); return QS(s.u * o.u + s.v * o.v * QS.d, s.u * o.v + s.v * o.u)
    __rmul__ = __mul__
    def inv(s):
        den = s.u * s.u - s.v * s.v * QS.d; return QS(s.u / den, -s.v / den)
    def __truediv__(s, o): o = o if isinstance(o, QS) else QS(o); return s * o.inv()
    def iszero(s): return s.u == 0 and s.v == 0
    def __repr__(s): return f'{s.u}' + (f' + {s.v} sqrt({QS.d})' if s.v else '')
def sqrt_rational(q):
    q = F(q)
    if q < 0: return None
    a, b = isqrt(q.numerator), isqrt(q.denominator)
    return F(a, b) if a * a == q.numerator and b * b == q.denominator else None
def make_a(d):
    """the element a = sqrt(d) in Q(sqrt d) (sets QS.d)"""
    QS.d = F(d); r = sqrt_rational(d)
    return QS(r, 0) if r is not None else QS(0, 1)
def solve_qs(cols, rhs, nrows_keys):
    """is rhs in span(cols)? cols: list of dicts key->QS; returns (consistent, rank)"""
    keys = nrows_keys; idx = {k: i for i, k in enumerate(keys)}
    M = [[QS(0)] * (len(cols) + 1) for _ in keys]
    for j, c in enumerate(cols):
        for k, v in c.items(): M[idx[k]][j] = M[idx[k]][j] + v
    for k, v in rhs.items(): M[idx[k]][len(cols)] = M[idx[k]][len(cols)] + v
    r = 0; ncol = len(cols)
    for c in range(ncol):
        piv = next((i for i in range(r, len(M)) if not M[i][c].iszero()), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]; inv = M[r][c].inv(); M[r] = [x * inv for x in M[r]]
        for i in range(len(M)):
            if i != r and not M[i][c].iszero():
                f = M[i][c]; M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        r += 1
    consistent = all(M[i][ncol].iszero() for i in range(r, len(M)))
    return consistent, r
# ======================= two-boson polynomial machinery (letters (comp, p), comp in {XC, FC}) =======================
def weight(m): return sum(p + 1 for (_, p) in m)
def monomials(w):
    letters = [(c, p) for c in (XC, FC) for p in range(w)]
    out = []
    def rec(start, rem, acc):
        if rem == 0: out.append(tuple(sorted(acc))); return
        for i in range(start, len(letters)):
            c, p = letters[i]
            if p + 1 <= rem: rec(i, rem - p - 1, acc + [letters[i]])
    rec(0, w, []); return out
def dpoly(Pl):
    out = {}
    for m, c in Pl.items():
        for i, (cc, p) in enumerate(m):
            nm = tuple(sorted(m[:i] + ((cc, p + 1),) + m[i + 1:])); out[nm] = out.get(nm, 0) + c
    return {k: v for k, v in out.items() if (v.iszero() is False if isinstance(v, QS) else v != 0)}
def residue_exp(G, alpha):
    """w-residue of G(z) exp(alpha.phi)(w) / exp(alpha.phi)(w), as a QS-polynomial in letters at w; alpha = {XC: QS, FC: QS}"""
    out = {}
    for m, c in G.items():
        n = len(m)
        for k in range(1, n + 1):
            for S in combinations(range(n), k):
                fac = QS(c) if not isinstance(c, QS) else c; P = 0
                for i in S:
                    cc, p = m[i]; fac = fac * alpha[cc] * ((-1) ** p * factorial(p)); P += p + 1
                T = P - 1
                rem = [m[i] for i in range(n) if i not in S]
                for sh in O._taylor(len(rem), T):
                    cf = fac
                    for kk in sh: cf = cf * F(1, factorial(kk))
                    nm = tuple(sorted((cc, p + kk) for (cc, p), kk in zip(rem, sh)))
                    out[nm] = out.get(nm, QS(0)) + cf
    return {k: v for k, v in out.items() if not v.iszero()}
def D_exp(Qm, alpha):
    """D(Q) = dQ + (alpha_x a_x + alpha_f a_f) Q for a single monomial Q (coefficient 1)"""
    out = {}
    for nm, v in dpoly({Qm: F(1)}).items(): out[nm] = out.get(nm, QS(0)) + QS(v)
    for cc in (XC, FC):
        nm = tuple(sorted(Qm + ((cc, 0),))); out[nm] = out.get(nm, QS(0)) + alpha[cc]
    return {k: v for k, v in out.items() if not v.iszero()}
def conserved(G, alpha):
    W = max(weight(m) for m in G)
    R = residue_exp(G, alpha)
    Qs = monomials(W - 2); cols = [D_exp(q, alpha) for q in Qs]
    keys = sorted(set(R) | {k for c in cols for k in c})
    ok, _ = solve_qs(cols, R, keys)
    return ok
def is_total_derivative(G):
    W = max(weight(m) for m in G); cols = [{k: QS(v) for k, v in dpoly({q: F(1)}).items()} for q in monomials(W - 1)]
    keys = sorted(set(G) | {k for c in cols for k in c}); ok, _ = solve_qs(cols, {k: QS(v) for k, v in G.items()}, keys); return ok
# ======================= Vir words =======================
def partitions(w, minpart):
    res = []
    def rec(rem, maxp, acc):
        if rem == 0: res.append(tuple(acc)); return
        for p_ in range(min(rem, maxp), minpart - 1, -1): rec(rem - p_, p_, acc + [p_])
    rec(w, w, []); return res
def words_pairs(Nval, w2, cache):
    T = {((), ((1, 1),)): F(1, 2)}
    def dT(k):
        f = T
        for _ in range(k): f = O.dfield_pairs(f)
        return f
    out = []
    for part in partitions(w2, 2):
        if ('p', Nval, part) not in cache:
            f = dT(part[-1] - 2)
            for n_ in reversed(part[:-1]): f = O.ope_pairs(dT(n_ - 2), f, Nval, (0,)).get(0, {})
            cache[('p', Nval, part)] = f
        out.append((part, cache[('p', Nval, part)]))
    return out
def words_phi(rho, w2, cache):
    T = {((FC, 0), (FC, 0)): F(1, 2), ((FC, 1),): F(rho)}
    def dT(k):
        f = T
        for _ in range(k): f = Wk.dfield(f)
        return f
    out = []
    for part in partitions(w2, 2):
        if ('f', rho, part) not in cache:
            f = dT(part[-1] - 2)
            for n_ in reversed(part[:-1]): f = Wk.ope(dT(n_ - 2), f, poles=(0,)).get(0, {})
            cache[('f', rho, part)] = f
        out.append((part, cache[('f', rho, part)]))
    return out
def pair_products(w):
    pairs = [(a, b) for a in range(1, w + 1) for b in range(a, w + 1) if a + b <= w]
    res = []
    def rec(start, rem, acc):
        if rem == 0: res.append(tuple(acc)); return
        for i in range(start, len(pairs)):
            a, b = pairs[i]
            if a + b <= rem: rec(i, rem - a - b, acc + [(a, b)])
    rec(0, w, []); return res
# ======================= exact rational solve with nullspace (pair coordinates) =======================
def solve_frac(cols, rhs):
    keys = sorted(set(rhs) | {k for c in cols for k in c}, key=str); idx = {k: i for i, k in enumerate(keys)}
    n = len(cols); M = [[F(0)] * (n + 1) for _ in keys]
    for j, c in enumerate(cols):
        for k, v in c.items(): M[idx[k]][j] += F(v)
    for k, v in rhs.items(): M[idx[k]][n] += F(v)
    piv_cols = []; r = 0
    for c in range(n):
        p_ = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if p_ is None: continue
        M[r], M[p_] = M[p_], M[r]; iv = 1 / M[r][c]; M[r] = [x * iv for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c]; M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        piv_cols.append(c); r += 1
    consistent = all(M[i][n] == 0 for i in range(r, len(M)))
    x = [F(0)] * n
    for i, c in enumerate(piv_cols): x[c] = M[i][n]
    free = [c for c in range(n) if c not in piv_cols]
    null = []
    for fc in free:
        v = [F(0)] * n; v[fc] = F(1)
        for i, c in enumerate(piv_cols): v[c] = -M[i][fc]
        null.append(v)
    return consistent, x, null
# ======================= record inputs =======================
def record_reps(tv, sol, weights):
    ROOT = os.path.join(HERE, '..', '..', '..', '..', 'code'); sys.path.insert(0, ROOT); cwd = os.getcwd(); os.chdir(ROOT)  # KIT PATCH: ROOT was an absolute path to a read-only clone of code/ (identical files)
    try:
        from invariant_engine import curve_point
        from mixed_engine import cyl_P4, gen_mixed_basis, residue_mixed_cur, d_mixed_mono, spin_mixed
        from sparse_linalg import genuine_kernel_sparse
        N_, s_ = curve_point(F(tv)); P4 = cyl_P4(sol, N_, s_); reps = {4: P4}
        for W in weights:
            if W == 4: continue
            g = genuine_kernel_sparse(P4, gen_mixed_basis(W), N_, residue_mixed_cur, gen_mixed_basis, d_mixed_mono, spin_mixed)
            assert len(g) == 1; reps[W] = g[0]
        return N_, s_, {W: {k: F(v) for k, v in c.items()} for W, c in reps.items()}
    finally: os.chdir(cwd)
# ======================= decomposition and translation =======================
def decompose_translate(rep, W, Nval, rho, cache):
    labels = []; cols = []
    for w2 in [0] + list(range(2, W + 1)):
        wds = [((), {((), ()): F(1)})] if w2 == 0 else words_pairs(Nval, w2, cache)
        for xpart in (partitions(W - w2, 1) if W - w2 > 0 else [()]):
            for part, f in wds:
                labels.append((xpart, part)); cols.append(O.mulx(xpart, f))
    nB = len(cols)
    for wy in range(0, W):
        wx = W - 1 - wy
        if wy == 1: continue
        for pp in ([()] if wy == 0 else pair_products(wy)):
            for xpart in (partitions(wx, 1) if wx > 0 else [()]):
                cols.append(O.dfield_pairs({(tuple(sorted(xpart)), tuple(pp)): F(1)}))
    ok, x, null = solve_frac(cols, rep)
    amb = [v for v in null if any(v[j] != 0 for j in range(nB))]
    def phi_image(coefs):
        G = {}
        for j in range(nB):
            if coefs[j] == 0: continue
            xpart, part = labels[j]
            word = {(): F(1)} if part == () else dict(words_phi(rho, sum(part), cache))[part]
            xl = tuple((XC, a - 1) for a in xpart)
            for m, c in word.items():
                nm = tuple(sorted(m + xl)); G[nm] = G.get(nm, 0) + coefs[j] * c
        return {k: v for k, v in G.items() if v != 0}
    G = phi_image(x)
    amb_td = [is_total_derivative(phi_image(v)) if phi_image(v) else True for v in amb]
    return ok, G, len(nB_ := [0] * nB), len(cols) - nB, len(null), len(amb), all(amb_td)
# ======================= run =======================
PAIRS = {1: [('pair', lambda p: (-p / 2, 1 - p * p / 4))], 2: [('pair', lambda p: (-p, 1 - p * p))],
         3: [('pair A', lambda p: (1 / p, 1 / (p * p))), ('pair B', lambda p: (-p, p * p))]}
FIBRES = {1: [F(8, 5), F(6, 5), F(40, 29)], 2: [F(4, 5), F(3, 5), F(8, 17)], 3: [F(1, 2), F(2, 3), F(3, 2)]}
OUT = {}; cache = {}
TRANSLATED = {}
for sol in (1, 2, 3):
    for fi, p in enumerate(FIBRES[sol]):
        tv = (2 + p * p) / (2 - p * p); rho = p / 2 - 1 / p
        t1 = time.time(); N_, s_, reps = record_reps(tv, sol, (4, 6, 8))
        assert N_ == 1 - 12 * rho * rho, 'c(rho) != N(t)'
        say(f'== Sol {sol}, p = {p}: t = {tv}, N = {N_} (= 1 - 12 rho^2, rho = {rho}), s = {s_}  (inputs {time.time()-t1:.0f} s)')
        for W in (4, 6, 8):
            t2 = time.time()
            ok, G, nB, nD, nnull, namb, ambtd = decompose_translate(reps[W], W, N_, rho, cache)
            result(f'VIR1 at Sol {sol} p = {p}, I_{W-1}: decomposition residual 0 (X x Vir {nB} + dV {nD}; nullspace {nnull}, ambiguous {namb}, ambiguity phi-images total derivatives {ambtd})', ok and ambtd)
            TRANSLATED[(sol, fi, W)] = G
            for name, fn in PAIRS[sol]:
                beta, d = fn(p); a = make_a(d)
                cp = conserved(G, {XC: a, FC: QS(beta)}) and conserved(G, {XC: -a, FC: QS(beta)})
                OUT[f'{sol}|{p}|{W}|{name}'] = cp
                result(f'Sol {sol} p = {p} I_{W-1} vs {name} (beta = {beta}, a^2 = {d}): V_+ and V_- conserved', cp)
            for bet in (p, -2 / p):
                QS.d = F(1)
                control(f'(c) Virasoro screening exp({bet} phi) conserved for I_{W-1} (Sol {sol}, p = {p})', conserved(G, {XC: QS(0), FC: QS(bet)}))
            say(f'      ({time.time()-t2:.0f} s)')
        # controls on I_5 at the first fibre
        if fi == 0:
            G5 = TRANSLATED[(sol, 0, 6)]
            for name, fn in PAIRS[sol]:
                beta, d = fn(p)
                a = make_a(d + F(1, 10)); bad = conserved(G5, {XC: a, FC: QS(beta)}) and conserved(G5, {XC: -a, FC: QS(beta)})
                control(f'(a) a^2 + 1/10 ({name}, Sol {sol}) FAILS on I_5', not bad)
            a = make_a(F(2, 7)); bad = conserved(G5, {XC: a, FC: QS(F(3, 11))}) and conserved(G5, {XC: -a, FC: QS(F(3, 11))})
            control(f'(d) generic (a^2, beta) = (2/7, 3/11) FAILS on I_5 (Sol {sol})', not bad)
            okm, Gm, *_ = decompose_translate(reps[6], 6, N_, -rho, cache)
            beta, d = PAIRS[sol][0][1](p); a = make_a(d)
            bad = conserved(Gm, {XC: a, FC: QS(beta)}) and conserved(Gm, {XC: -a, FC: QS(beta)})
            control(f'(b) rho -> -rho translation FAILS for {PAIRS[sol][0][0]} on I_5 (Sol {sol})', not bad)
            QS.d = F(1)
            san = {k: v for k, v in dpoly({((XC, 0), (XC, 0), (FC, 2)): F(1), ((FC, 0),) * 3: F(5, 7)}).items()}
            beta, d = PAIRS[sol][0][1](p); a = make_a(d)
            control(f'(f) sanity: a total derivative passes V_+- (Sol {sol})', conserved(san, {XC: a, FC: QS(beta)}) and conserved(san, {XC: -a, FC: QS(beta)}))
json.dump({k: v for k, v in OUT.items()}, open(os.path.join(HERE, 'vir2r_conservation.json'), 'w'), indent=1)
json.dump({f'{s}|{fi}|{W}': {str(m): str(c) for m, c in G.items()} for (s, fi, W), G in TRANSLATED.items()}, open(os.path.join(HERE, 'vir2r_translated.json'), 'w'), indent=1)
say(f'conservation phase done ({time.time() - T0:.0f} s)')
# ======================= commutant (first fibre per solution) =======================
def commutant_dim(W, exps):
    """dim of {G weight W : R_E(G) = D_E(Q_E) for all E} modulo d(weight W-1); exact over Q (all exps have rational parameters here)"""
    Gm = monomials(W); Qm = monomials(W - 2); nG, nQ = len(Gm), len(Qm)
    rows = {}
    ncol = nG + nQ * len(exps)
    for e, alpha in enumerate(exps):
        for j, g in enumerate(Gm):
            for k, v in residue_exp({g: F(1)}, alpha).items():
                rows.setdefault((e, k), {})[j] = rows.get((e, k), {}).get(j, F(0)) + (v.u if v.v == 0 else None)
        for j, q in enumerate(Qm):
            for k, v in D_exp(q, alpha).items():
                c = nG + e * nQ + j; rows.setdefault((e, k), {})[c] = rows.get((e, k), {}).get(c, F(0)) - v.u
    # nullspace via elimination (Fractions)
    R = [dict(r) for r in rows.values()]
    piv = {}; order = []
    for r in R:
        r = {k: v for k, v in r.items() if v != 0}
        while r:
            k0 = min(r)
            if k0 in piv:
                f = r[k0]; pr = piv[k0]
                for k, v in pr.items(): r[k] = r.get(k, F(0)) - f * v
                r = {k: v for k, v in r.items() if v != 0}
            else:
                iv = 1 / r[k0]; piv[k0] = {k: v * iv for k, v in r.items()}; order.append(k0); break
    # back-substitute to reduced form
    for k0 in sorted(piv, reverse=True):
        pr = piv[k0]
        for k1 in list(piv):
            if k1 != k0 and k0 in piv[k1]:
                f = piv[k1][k0]; piv[k1] = {k: piv[k1].get(k, F(0)) - f * pr.get(k, F(0)) for k in set(piv[k1]) | set(pr)}
                piv[k1] = {k: v for k, v in piv[k1].items() if v != 0}
    free = [c for c in range(ncol) if c not in piv]
    basisG = []
    for fc in free:
        v = {fc: F(1)}
        for k0, pr in piv.items():
            if fc in pr: v[k0] = -pr[fc]
        g = {Gm[c]: val for c, val in v.items() if c < nG and val != 0}
        basisG.append(g)
    # rank of the G-projections
    keys = {m: i for i, m in enumerate(Gm)}
    rank_rows = []
    pv = {}
    rk = 0
    for g in basisG:
        r = {keys[m]: v for m, v in g.items()}
        while r:
            k0 = min(r)
            if k0 in pv:
                f = r[k0]
                for k, v in pv[k0].items(): r[k] = r.get(k, F(0)) - f * v
                r = {k: v for k, v in r.items() if v != 0}
            else:
                iv = 1 / r[k0]; pv[k0] = {k: v * iv for k, v in r.items()}; rk += 1; break
    nd = len(monomials(W - 1))
    return rk - nd, rk, nd
if not NOCOMM:
    for sol in (1, 2, 3):
        p = FIBRES[sol][0]; beta, d = PAIRS[sol][0][1](p); a0 = sqrt_rational(d); assert a0 is not None
        QS.d = F(1)
        exps = [{XC: QS(0), FC: QS(p)}, {XC: QS(a0), FC: QS(beta)}, {XC: QS(-a0), FC: QS(beta)}]
        for W in (4, 6, 8):
            t3 = time.time(); dimc, rk, nd = commutant_dim(W, exps)
            G = TRANSLATED[(sol, 0, W)]
            notd = not is_total_derivative(G)
            result(f'COMMUTANT Sol {sol} p = {p} weight {W}: dim mod d = {dimc} (proj rank {rk} - d-image {nd}); == 1 and spanned by I_{W-1} (in commutant by conservation; not a total derivative: {notd})',
                   dimc == 1 and notd and OUT.get(f'{sol}|{p}|{W}|{PAIRS[sol][0][0]}', False))
            say(f'      ({time.time()-t3:.0f} s)')
say(f'done ({time.time() - T0:.0f} s)')
if FAIL: say(f'CONTROLS FAILED {FAIL} -> exit 1'); sys.exit(1)
if RES: say(f'SEALED TESTS FAILED {RES} -> exit 2 (a result)'); sys.exit(2)
say('controls passed; sealed prediction holds -> exit 0'); sys.exit(0)
