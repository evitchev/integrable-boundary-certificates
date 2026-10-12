"""WUW1 -- object: at t = 5/3 (p^2 = 1/2), bosons U, W = (phi +- X)/sqrt2.  W_U := Res e^(2U)(z) e^(-2U)(w) (c_U = -2 Coulomb gas, alpha = 2, -1);
G_6 := W_U W_W.  Tests (sealed SEAL_WUW1.md 425bf68b...): W_U properties; I_3 = weight-4 commutant of {e^((U+W)/2), e^(2U), e^(2W)}; K_6 = ker ad I_3 at
weight 6 mod d; G_6 in K_6?  eigen(W_U) ~ P_U(P_U^2 - 1/4)?  eigen(G_6) ~ E_1?  Gates vs cyl(3,4), vev_sol3_w6, vev_sol2_w6.  Exact over Q."""
import os, sys, json, hashlib, time
from fractions import Fraction as F
from math import factorial
from itertools import combinations
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
assert hashlib.sha256(open(os.path.join(HERE, 'SEAL_WUW1.md'), 'rb').read()).hexdigest().startswith('425bf68b')
ORIG = {l.split()[1]: l.split()[0] for l in open(os.path.join(HERE, 'inputs', 'ORIGIN_SHA256'))}
assert all(hashlib.sha256(open(os.path.join(HERE, 'inputs', f_), 'rb').read()).hexdigest() == h for f_, h in ORIG.items())
import sympy as sp
import wick as Wk
LOG = open(os.path.join(HERE, 'run_wuw1.log'), 'w'); T0 = time.time()
def say(*a):
    s_ = ' '.join(str(x) for x in a); print(s_, flush=True); LOG.write(s_ + '\n'); LOG.flush()
FAIL = []; RES = {}
def control(n_, ok):
    say(f'CONTROL {n_}: {ok}')
    if not ok: FAIL.append(n_)
def result(n_, ok):
    say(f'SEALED TEST {n_}: {ok}'); RES[n_] = bool(ok)
U, Wc = 1, 2
# ---------------- basics
def weight(m): return sum(q + 1 for (_, q) in m)
def monomials(w):
    letters = [(c, q) for c in (U, Wc) for q in range(w)]
    out = []
    def rec(start, rem, acc):
        if rem == 0: out.append(tuple(sorted(acc))); return
        for i in range(start, len(letters)):
            c, q = letters[i]
            if q + 1 <= rem: rec(i, rem - q - 1, acc + [(c, q)])
    rec(0, w, []); return sorted(set(out))
def dpoly(G): return Wk.dfield(G)
def add(A, B, s=1): return Wk.fadd(A, B, s)
def scale(A, s): return Wk.fscale(A, s)
def mul(A, B):
    out = {}
    for m1, c1 in A.items():
        for m2, c2 in B.items():
            m = tuple(sorted(m1 + m2)); out[m] = out.get(m, 0) + c1 * c2
    return {m: c for m, c in out.items() if c}
# ---------------- exact linear algebra over Q
def rref(rows, ncol):
    M = [list(r) for r in rows]; piv = []; r = 0
    for c in range(ncol):
        p = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]; inv = 1 / M[r][c]; M[r] = [x * inv for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c]; M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        piv.append(c); r += 1
        if r == len(M): break
    return M[:r], piv
def nullspace(cols_as_vectors, nrows):
    """columns c_j (dict row->val); returns basis of {x : sum x_j c_j = 0}"""
    ncol = len(cols_as_vectors)
    rows = [[F(0)] * ncol for _ in range(nrows)]
    for j, col in enumerate(cols_as_vectors):
        for i, v in col.items(): rows[i][j] += v
    R, piv = rref(rows, ncol); free = [j for j in range(ncol) if j not in piv]; basis = []
    for f in free:
        x = [F(0)] * ncol; x[f] = F(1)
        for ri, pc in enumerate(piv): x[pc] = -R[ri][f]
        basis.append(x)
    return basis
def rank_of(vectors, keys):
    idx = {k: i for i, k in enumerate(keys)}
    rows = [[F(0)] * len(keys) for _ in vectors]
    for r, v in enumerate(vectors):
        for k, c in v.items(): rows[r][idx[k]] += c
    R, piv = rref(rows, len(keys)); return len(piv)
def mod_d_dim(fields, W):
    """dim of span(fields) modulo d(weight W-1 monomials)"""
    dims = [dpoly({m: F(1)}) for m in monomials(W - 1)]
    keys = sorted({k for f in fields + dims for k in f})
    return rank_of(fields + dims, keys) - rank_of(dims, keys)
def in_span_mod_d(target, fields, W):
    dims = [dpoly({m: F(1)}) for m in monomials(W - 1)]; keys = sorted({k for f in fields + dims + [target] for k in f})
    return rank_of(fields + dims + [target], keys) == rank_of(fields + dims, keys)
# ---------------- exponentials: zero-mode test (residue mod D) and field-level screening
def residue_exp(G, alpha):
    """Res_z G(z) e^(alpha.Phi)(w) / e^(alpha.Phi)(w) as a polynomial in letters at w"""
    out = {}
    for m, c in G.items():
        n = len(m)
        for k in range(1, n + 1):
            for S in combinations(range(n), k):
                fac = F(c); P = 0
                for i in S:
                    cc, q = m[i]; fac *= alpha[cc] * ((-1) ** q * factorial(q)); P += q + 1
                rem = [m[i] for i in range(n) if i not in S]
                for sh, cf in Wk._taylor(rem, P - 1):
                    nm = tuple(sorted(sh)); out[nm] = out.get(nm, 0) + fac * cf
    return {k: v for k, v in out.items() if v}
def D_exp(Qm, alpha):
    out = dict(dpoly({Qm: F(1)}))
    for cc in (U, Wc):
        if alpha[cc]: nm = tuple(sorted(Qm + ((cc, 0),))); out[nm] = out.get(nm, 0) + alpha[cc]
    return {k: v for k, v in out.items() if v}
def schur(alpha, comp, n):
    """coefficient of x^n in exp(sum_k alpha d^k Phi_comp x^k / k!) as a field (letters (comp, k-1))"""
    ser = {0: {(): F(1)}}
    gen = {k: {((comp, k - 1),): F(alpha) / factorial(k)} for k in range(1, n + 1)}
    # exp via recursion: n e_n = sum_k k g_k e_(n-k)
    e = [{(): F(1)}]
    for m in range(1, n + 1):
        acc = {}
        for k in range(1, m + 1): acc = add(acc, scale(mul(gen[k], e[m - k]), F(k)))
        e.append(scale(acc, F(1, m)))
    return e[n]
def field_screen(G, alpha_c, comp):
    """Res_z e^(alpha U_comp)(z) G(w) / e^(alpha U_comp)(w): field-level action of the screening on G"""
    out = {}
    for m, c in G.items():
        n = len(m); idxs = [i for i in range(n) if m[i][0] == comp]
        for k in range(0, len(idxs) + 1):
            for S in combinations(idxs, k):
                fac = F(c); P = 0
                for i in S:
                    q = m[i][1]; fac *= -alpha_c * factorial(q); P += q + 1
                if P == 0: continue
                rem = tuple(m[i] for i in range(n) if i not in S)
                sc = schur(alpha_c, comp, P - 1)
                out = add(out, scale(mul(sc, {tuple(sorted(rem)): F(1)}), fac))
    return out
def commutant(W, alphas):
    """basis (fields) of {G weight W : Res(G, alpha) = D_alpha(Q_alpha) for all alpha} (not yet mod d)"""
    gm = monomials(W); qm = monomials(W - 2)
    cols = []; rowkeys = {}
    def key(a_i, mono):
        k_ = (a_i, mono)
        if k_ not in rowkeys: rowkeys[k_] = len(rowkeys)
        return rowkeys[k_]
    for g in gm:
        col = {}
        for ai, al in enumerate(alphas):
            for mono, v in residue_exp({g: F(1)}, al).items(): col[key(ai, mono)] = col.get(key(ai, mono), 0) + v
        cols.append(col)
    for ai, al in enumerate(alphas):
        for q in qm:
            col = {}
            for mono, v in D_exp(q, al).items(): col[key(ai, mono)] = col.get(key(ai, mono), 0) - v
            cols.append(col)
    ns = nullspace(cols, len(rowkeys))
    return [{gm[j]: x[j] for j in range(len(gm)) if x[j] != 0} for x in ns if any(x[j] != 0 for j in range(len(gm)))]
def reduce_mod_d(fields, W):
    """independent representatives mod d"""
    reps = []
    for f in fields:
        if mod_d_dim(reps + [f], W) > mod_d_dim(reps, W) if reps else mod_d_dim([f], W) > 0: reps.append(f)
    return reps
def ad_kernel(I, W):
    """{G weight W : Res_1 I(z) G(w) in d(...)} mod d"""
    gm = monomials(W); wI = max(weight(m) for m in I); wr = wI + W - 1
    dm = monomials(wr - 1); dcols = [dpoly({m: F(1)}) for m in dm]
    rowkeys = {}
    def key(mono):
        if mono not in rowkeys: rowkeys[mono] = len(rowkeys)
        return rowkeys[mono]
    cols = []
    for g in gm:
        r1 = Wk.ope(I, {g: F(1)}, poles={1}).get(1, {})
        cols.append({key(k): v for k, v in r1.items()})
    for dc in dcols: cols.append({key(k): -v for k, v in dc.items()})
    ns = nullspace(cols, len(rowkeys))
    fields = [{gm[j]: x[j] for j in range(len(gm)) if x[j] != 0} for x in ns if any(x[j] != 0 for j in range(len(gm)))]
    return reduce_mod_d(fields, W)
# ---------------- cylinder recipe (VIR3-R, own)
eps = sp.Symbol('eps')
Rser = sp.series(1 / (4 * sp.sinh(eps / 2)**2) - 1 / eps**2, eps, 0, 26).removeO()
Rd = [sp.Rational(sp.diff(Rser, eps, k).subs(eps, 0)) for k in range(24)]
def Cab(a, b): return (-1)**b * Rd[a + b]
PU, PW, PX, PI = sp.symbols('P_U P_W P_X pi')
def pairings(L):
    if not L: yield [], []; return
    first, rest = L[0], L[1:]
    for pr, un in pairings(rest): yield pr, [first] + un
    for i, o in enumerate(rest):
        if o[0] == first[0]:
            for pr, un in pairings(rest[:i] + rest[i + 1:]): yield [(first, o)] + pr, un
def eigen(G, internal=True):
    tot = 0
    for m, c in G.items():
        for pr, un in pairings(list(m)):
            if any(q >= 1 for (_, q) in un): continue
            val = sp.Rational(c.numerator, c.denominator)
            for (a, b) in pr: val *= Cab(a[1], b[1])
            for (cc, q) in un: val *= (PU if cc == U else PW)
            tot += val
    tot = sp.expand(tot)
    return sp.expand(tot.subs({PU: (PX + PI) / sp.sqrt(2), PW: (PI - PX) / sp.sqrt(2)})) if internal else tot
# ---------------- certified side at t = 5/3
tv = sp.Rational(5, 3); p2 = sp.Rational(1, 2); rho2 = (p2 - 2)**2 / (4 * p2)
Pa, Qa, tS = sp.symbols('P Q t')
ROOT = os.path.join(HERE, '..', '..', '..', '..', 'code'); sys.path.insert(0, ROOT); _cwd = os.getcwd(); os.chdir(ROOT)  # KIT PATCH: repo code/ from __file__ (was ~/opus5-w8-clone/code; cyl_first_order.py, wkb.py, first_order/vev_w4_points.json byte-identical there)
from cyl_first_order import cyl, X as Xs, Y as Ys, t as tc
os.chdir(_cwd)
def to_internal(e_art):
    return sp.expand(e_art.subs({Pa: sp.I * PX}).subs({Qa: sp.sqrt(rho2 - PI**2)}))
def table(sol, Wt):
    d = json.load(open(os.path.join(HERE, 'inputs', f'vev_sol{sol}_w{Wt}.json')))
    return to_internal(sp.expand(sum(sp.cancel(sp.sympify(v, locals={'t': tS}).subs(tS, tv)) * Pa**int(k_.split()[0][2:]) * Qa**int(k_.split()[1][2:]) for k_, v in d['coefficients'].items())))
CERT3 = {3: to_internal(sp.expand(cyl(3, 4).subs(tc, tv).subs({Xs: Pa**2, Ys: Qa**2}))), 5: table(3, 6)}
CERT2_5 = table(2, 6)
def proportional(a, b):
    if a == 0 or b == 0: return False
    pa, pb = sp.Poly(a, PX, PI).as_dict(), sp.Poly(b, PX, PI).as_dict(); mons = set(pa) | set(pb)
    ref = next(iter(pa)); return all(sp.simplify(pa.get(m, 0) * pb.get(ref, 0) - pb.get(m, 0) * pa[ref]) == 0 for m in mons) and pb.get(ref, 0) != 0
def in_eig_span(target, vecs):
    syms = sp.symbols(f'c0:{len(vecs)}'); expr = sp.expand(target - sum(s_ * v for s_, v in zip(syms, vecs)))
    eqs = sp.Poly(expr, PX, PI).coeffs(); sol = sp.solve(eqs, syms, dict=True); return (sol[0] if sol else None)
# ================= (1) the weight-3 objects
TU = {((U, 0), (U, 0)): F(1, 2), ((U, 1),): F(1, 2)}
def Wcur(a, comp): return schur(F(a), comp, 3)
WU, WW = Wcur(2, U), Wcur(2, Wc)
say(f'(1) W_U = Res e^(2U)(z) e^(-2U)(w) = {WU}')
o = Wk.ope(TU, WU)
prim = (all(not o.get(r) for r in o if r >= 3) and o.get(2, {}) == scale(WU, F(3)) and o.get(1, {}) == dpoly(WU))
control('(1)(i) W_U is T_U-primary of weight 3 (c_U = -2)', prim)
say(f'   c(T_U) check: T_U T_U 4th pole = {Wk.ope(TU, TU).get(4)}  (c/2)')
control('(1)(ii) field level: Res e^(-U)(z) W_U(w) = 0 (short screening)', not field_screen(WU, F(-1), U))
say(f'   (reported) field level Res e^(2U)(z) W_U(w) = {field_screen(WU, F(2), U)}')
dT = dpoly(TU); notdesc = rank_of([dT, WU], sorted(set(dT) | set(WU))) == 2
control('(1)(iii) W_U not a T_U descendant (not in span{dT_U})', notdesc)
eWU = sp.factor(eigen(WU, internal=False)); say(f'   eigen(W_U) = {eWU}')
result('P1 eigen(W_U) proportional to P_U (P_U^2 - 1/4)', sp.simplify(sp.cancel(eWU / (PU * (PU**2 - sp.Rational(1, 4))))).free_symbols == set())
# ================= (2) I_3, K_6
SCR = [{U: F(1, 2), Wc: F(1, 2)}, {U: F(2), Wc: F(0)}, {U: F(0), Wc: F(2)}]
c4 = reduce_mod_d(commutant(4, SCR), 4); say(f'(2) weight-4 commutant of {{e^((U+W)/2), e^(2U), e^(2W)}} mod d: dim {len(c4)}')
I3 = c4[0]; e3 = eigen(I3)
control('GATE I_3: dim 1 and eigen(I_3) proportional to the record cyl(3,4) at t = 5/3', len(c4) == 1 and proportional(e3, CERT3[3]))
c6 = reduce_mod_d(commutant(6, SCR), 6); say(f'   weight-6 screening commutant mod d: dim {len(c6)}')
T1 = time.time(); K6 = ad_kernel(I3, 6); say(f'   K_6 = ker ad I_3 at weight 6 mod d: dim {len(K6)} ({time.time()-T1:.0f}s)')
eK = [eigen(k_) for k_ in K6]
control('GATE K_6 contains the weight-6 screening commutant (I_5)', all(in_span_mod_d(c, K6, 6) for c in c6))
control('GATE vev_sol3_w6 (t = 5/3) in span eigen(K_6)', in_eig_span(CERT3[5], eK) is not None)
s2 = in_eig_span(CERT2_5, eK); control('GATE vev_sol2_w6 (t = 5/3) in span eigen(K_6)', s2 is not None)
control('eigen of a total derivative vanishes (recipe sanity)', eigen(dpoly(WU)) == 0 and eigen(dpoly(mul(TU, WU))) == 0)
G6 = mul(WU, WW)
inK = in_span_mod_d(G6, K6, 6); notI5 = not in_span_mod_d(G6, c6, 6)
result('P2 G_6 = W_U W_W in K_6 mod d, and not in span{I_5} + d', inK and notI5)
say(f'   (P2 parts: in K_6 {inK}; not in I_5 + d {notI5})')
eG = eigen(G6); E1 = sp.expand(sp.prod([(PX + PI - j * sp.sqrt(p2)) * (PX - PI - j * sp.sqrt(p2)) for j in (-1, 0, 1)]))
say(f'   eigen(G_6) = {sp.factor(eG)}')
say(f'   E_1 = {sp.factor(E1)}')
result('P3 eigen(G_6) proportional to E_1', proportional(E1, eG))
anti = sp.expand(eG.subs({PX: PI, PI: PX}, simultaneous=True) + eG) == 0
control('antisymmetry: eigen(G_6) odd under P_X <-> pi', anti)
say(f"   Fable's relation at 5/3 with this E: Sol 2 in span(Sol 3, eigen(G_6))? {in_eig_span(CERT2_5, [CERT3[5], eG])}")
for nm, al in (('pair B e^(-U)', {U: F(-1), Wc: F(0)}), ('pair B e^(-W)', {U: F(0), Wc: F(-1)}), ('Virasoro e^((U+W)/2)', SCR[0]), ('pair A e^(2U)', SCR[1]), ('Virasoro e^(-2(U+W))', {U: F(-2), Wc: F(-2)})):
    gm = monomials(4); R = residue_exp(G6, al); cols = [D_exp(q, al) for q in gm]; keys = sorted(set(R) | {k for c in cols for k in c})
    ok = rank_of(cols + [R], keys) == rank_of(cols, keys)
    say(f'   P4 (reported): G_6 zero mode commutes with {nm}: {ok}')
# ================= (3) controls
Gs = add(mul(WU, WU), mul(WW, WW)); eS = eigen(Gs)
control('symmetric W_U^2 + W_W^2: eigenvalue NOT antisymmetric', sp.expand(eS.subs({PX: PI, PI: PX}, simultaneous=True) + eS) != 0)
say(f'   (reported) W_U^2 + W_W^2 in K_6 mod d: {in_span_mod_d(Gs, K6, 6)}; eigen symmetric: {sp.expand(eS.subs({PX: PI, PI: PX}, simultaneous=True) - eS) == 0}')
Gt = add(G6, {((U, 0),) * 6: F(1, 100)})
control('TAMPER 1 G_6 + (dU)^6/100 leaves K_6 and breaks P3', (not in_span_mod_d(Gt, K6, 6)) and not proportional(E1, eigen(Gt)))
WUt = Wcur(F(21, 10), U); WWt = Wcur(F(21, 10), Wc)
control('TAMPER 2 exponent 21/10: breaks P1 and P3', sp.simplify(sp.cancel(eigen(WUt, internal=False) / (PU * (PU**2 - sp.Rational(1, 4))))).free_symbols != set() and not proportional(E1, eigen(mul(WUt, WWt))))
E1t = sp.expand(sp.prod([(PX + PI - j * sp.sqrt(p2 + sp.Rational(1, 1000))) * (PX - PI - j * sp.sqrt(p2 + sp.Rational(1, 1000))) for j in (-1, 0, 1)]))
control('SCORING TAMPER E_1 with p^2 + 1/1000 breaks P3', not proportional(E1t, eG))
json.dump({'W_U': {str(k): str(v) for k, v in WU.items()}, 'I3': {str(k): str(v) for k, v in I3.items()}, 'K6': [{str(k): str(v) for k, v in f.items()} for f in K6],
           'eigen_G6': str(eG), 'E1': str(E1), 'eigen_K6': [str(e) for e in eK], 'results': RES}, open(os.path.join(HERE, 'wuw1.json'), 'w'), indent=1)
say(f'done ({time.time()-T0:.0f}s)')
if FAIL: say(f'CONTROLS FAILED {FAIL} -> exit 1'); sys.exit(1)
if not all(RES.values()): say('a sealed prediction failed -> exit 2 (a result)'); sys.exit(2)
sys.exit(0)
