"""Library of the q = 8 density-level test. On top of the byte-copied WUW2 engine kwick.py/wick.py (in code/vendor/).
(i) the kappa-weighted cylinder vacuum rule exactly as wuw2.py's eigen_raw; (ii) the record dictionary; (iii) the Euler operators
E_V, E_Z (P in dF iff both vanish, weight >= 1); (iv) the complement basis C_w of dF_{w-1} in F_w; (v) exact linear algebra with
python-flint (fmpq_mat over Q; nmod_mat mod p for cross-checks); (vi) record readers; (vii) Virasoro vacuum-module fields."""
import os, sys, json, hashlib, time
from fractions import Fraction as F
from math import factorial
HERE = os.path.dirname(os.path.abspath(__file__)); WS = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HERE, 'vendor'))
import kwick as K
import wick as Wk
from kwick import A, B
import sympy as sp
import flint
REPO = os.path.abspath(os.path.join(WS, '..', '..', '..', '..'))   # repository root: this folder is results/lab/jumping_loci_records/q8_density_test/
def say(*a): print(' '.join(str(x) for x in a), flush=True)
def self_hash(path): return hashlib.sha256(open(path, 'rb').read()).hexdigest()
def p2(n):
    pp = [1] + [0] * n
    for k in range(1, n + 1):
        for j in range(k, n + 1): pp[j] += pp[j - k]
    return sum(pp[i] * pp[n - i] for i in range(n + 1))
# ---------------- (i) cylinder vacuum rule (wuw2.py eigen_raw, verbatim in substance)
_eps = sp.Symbol('eps')
_Rser = sp.series(1 / (4 * sp.sinh(_eps / 2)**2) - 1 / _eps**2, _eps, 0, 34).removeO()
RD = [sp.Rational(sp.diff(_Rser, _eps, k).subs(_eps, 0)) for k in range(32)]
PA, PB, PX, PI = sp.symbols('P_A P_B P_X pi')
def pairings(L):
    if not L: yield [], []; return
    first, rest = L[0], L[1:]
    for pr, un in pairings(rest): yield pr, [first] + un
    for i, o in enumerate(rest):
        if o[0] == first[0]:
            for pr, un in pairings(rest[:i] + rest[i + 1:]): yield [(first, o)] + pr, un
def eigen_raw(G):
    tot = 0; kap = sp.Rational(K.KAPPA.numerator, K.KAPPA.denominator)
    for m, c in G.items():
        for pr, un in pairings(list(m)):
            if any(q >= 1 for (_, q) in un): continue
            val = sp.Rational(c.numerator, c.denominator)
            for (a, b) in pr: val *= kap * (-1)**b[1] * RD[a[1] + b[1]]
            for (cc, q) in un: val *= (PA if cc == A else PB)
            tot += val
    return sp.expand(tot)
def internal(e, kappa):
    """P_A -> P_V, P_B -> P_Z in the record's (P_X, pi): V = U/s, s^2 = 1/kappa, P_U = (P_X + pi)/sqrt2 -> P_V = sqrt(kappa/2)(P_X + pi)."""
    r = sp.sqrt(sp.Rational(kappa) / 2)
    return sp.expand(e.subs({PA: r * (PX + PI), PB: r * (PI - PX)}, simultaneous=True))
def proportional(a, b, gens):
    if a == 0 or b == 0: return False
    pa, pb = sp.Poly(a, *gens).as_dict(), sp.Poly(b, *gens).as_dict(); ref = next(iter(pa))
    return pb.get(ref, 0) != 0 and all(sp.simplify(pa.get(m, 0) * pb[ref] - pb.get(m, 0) * pa[ref]) == 0 for m in set(pa) | set(pb))
def swap_sym(e): return sp.expand(e.subs({PX: PI, PI: PX}, simultaneous=True))
# ---------------- (ii) record dictionary and readers
def rho2_of_t(tv): return sp.Rational(2) / (sp.Rational(tv)**2 - 1)
def rec_expr(terms, tv):
    """sum_{a,b} v_ab X_rec^a Y_rec^b -> polynomial in (P_X, pi): X_rec = -P_X^2, Y_rec = rho^2 - pi^2."""
    r2 = rho2_of_t(tv); e = 0
    for (a_, b_), v in terms.items(): e += v * (-PX**2)**a_ * (r2 - PI**2)**b_
    return sp.expand(e)
def read_cyl34(sol, tstr):
    d = json.load(open(os.path.join(REPO, 'results/lab/first_order/vev_w4_points.json')))[str(sol)][tstr]
    return rec_expr({tuple(int(z) for z in k.split(',')): sp.Rational(v) for k, v in d.items()}, sp.Rational(tstr))
def read_vev_profile(path, tv):
    """anchor-type JSON, key 'vev': {'a,b': rational function of t} (as z5_gate_K6_K7.py)."""
    t = sp.Symbol('t'); tab = json.load(open(path))['vev']; terms = {}
    for key, val in tab.items():
        a_, b_ = [int(x) for x in key.split(',')]
        f = sp.sympify(val, locals={'t': t}, rational=True); num, den = sp.fraction(sp.together(f))
        assert den.subs(t, tv) != 0
        terms[(a_, b_)] = sp.Rational(num.subs(t, tv)) / sp.Rational(den.subs(t, tv))
    return rec_expr(terms, tv)
def read_vev_coeffs(path, tv):
    """vev_sol*_w*.json, key 'coefficients': {'P^i Q^j': rational function of t} in article momenta (as k10_full.py)."""
    t = sp.Symbol('t'); d = json.load(open(path)); Pa, Qa = sp.symbols('P Q'); e = 0
    for k_, v in d['coefficients'].items():
        i, j = [int(z) for z in k_.replace('P^', '').replace('Q^', '').split()]
        e += sp.cancel(sp.sympify(v, locals={'t': t}).subs(t, tv)) * Pa**i * Qa**j
    r2 = rho2_of_t(tv)
    return sp.expand(sp.expand(e).subs({Pa: sp.I * PX}).subs({Qa: sp.sqrt(r2 - PI**2)}))
# ---------------- (iii) Euler operators
_ECACHE = {}
def euler_mono(m, comp):
    key = (m, comp)
    if key in _ECACHE: return _ECACHE[key]
    out = {}; seen = set()
    for idx, (cc, q) in enumerate(m):
        if cc != comp or q in seen: continue
        seen.add(q); e = sum(1 for (c2, q2) in m if c2 == comp and q2 == q)
        rest = list(m[:idx] + m[idx + 1:]); fac = e * (-1)**q * factorial(q)
        for sh, cf in Wk._taylor(rest, q):
            mono = tuple(sorted(sh)); v = out.get(mono, 0) + fac * cf
            if v: out[mono] = v
            else: out.pop(mono, None)
    out = {mm: (int(v) if F(v).denominator == 1 else F(v)) for mm, v in out.items()}
    _ECACHE[key] = out; return out
def euler(P, comp):
    out = {}
    for m, c in P.items():
        for mm, v in euler_mono(m, comp).items():
            x = out.get(mm, 0) + c * v
            if x: out[mm] = x
            else: out.pop(mm, None)
    return out
def euler_vec(P):
    v = {}
    for comp in (A, B):
        for mm, x in euler(P, comp).items(): v[(comp, mm)] = x
    return v
def is_exact(P): return not euler(P, A) and not euler(P, B)
# ---------------- (iv) complement basis
def complement_basis(w):
    out = []
    for m in K.monomials(w):
        zs = [q for (c, q) in m if c == B]; vs = [q for (c, q) in m if c == A]
        if not zs:
            top = max(vs)
            if vs.count(top) >= 2: out.append(m)
        else:
            top = max(zs)
            if zs.count(top) >= 2 or (len(zs) == 1 and top == 0): out.append(m)
    return out
def check_complement(w):
    e = K.d_ech(w); n0 = len(e.rows); C = complement_basis(w); new = sum(1 for m in C if e.insert({m: F(1)}))
    return n0, len(C), new, len(e.rows), p2(w), p2(w - 1)
# ---------------- (v) exact linear algebra (python-flint)
def _fq(x):
    x = F(x); return flint.fmpq(x.numerator, x.denominator)
def build_mat(vectors, keys=None):
    if keys is None: keys = sorted({k for v in vectors for k in v})
    idx = {k: i for i, k in enumerate(keys)}; n = len(vectors)
    ent = [0] * (len(keys) * n)
    for j, v in enumerate(vectors):
        for k, x in v.items(): ent[idx[k] * n + j] = _fq(x)
    return flint.fmpq_mat(len(keys), n, ent), keys
def rank_Q(vectors):
    M, _ = build_mat(vectors); return M.rank()
def nullspace_Q(vectors):
    """basis of {x in Q^n : sum_j x_j v_j = 0}, from the rref over Q (lists of Fractions)."""
    M, keys = build_mat(vectors); n = len(vectors)
    if M.nrows() == 0: return [[F(1) if i == j else F(0) for i in range(n)] for j in range(n)], 0
    R, rk = M.rref(); piv = []
    for i in range(rk):
        j = next(j for j in range(n) if R[i, j] != 0); piv.append(j)
    free = [j for j in range(n) if j not in piv]; basis = []
    for f_ in free:
        x = [F(0)] * n; x[f_] = F(1)
        for i, pc in enumerate(piv):
            r = R[i, f_]; x[pc] = -F(int(r.p), int(r.q))
        basis.append(x)
    return basis, rk
def rank_mod(vectors, p):
    keys = sorted({k for v in vectors for k in v}); idx = {k: i for i, k in enumerate(keys)}; n = len(vectors)
    ent = [0] * (len(keys) * n)
    for j, v in enumerate(vectors):
        for k, x in v.items():
            x = F(x); assert x.denominator % p != 0
            ent[idx[k] * n + j] = (x.numerator % p) * pow(x.denominator, -1, p) % p
    return flint.nmod_mat(len(keys), n, ent, p).rank()
def consistent_Q(vectors, rhs):
    """is rhs in span(vectors)? by ranks over Q"""
    r0 = rank_Q(vectors); r1 = rank_Q(vectors + [rhs]); return r0 == r1, r0, r1
# ---------------- (vi) fields
def poles_vec(T, G):
    o = K.ope(T, G); v = {}
    for r, f in o.items():
        if r >= 3:
            for m, c in f.items(): v[(r, m)] = v.get((r, m), 0) + c
    return {k: x for k, x in v.items() if x}
def is_primary(T, G, h):
    o = K.ope(T, G)
    return all(not o.get(r) for r in o if r >= 3) and o.get(2, {}) == K.scale(G, F(h)) and o.get(1, {}) == K.dpoly(G)
def vir_vacuum_fields(T, w):
    """fields of L_{-n1}...L_{-nk}|0>, n1 >= ... >= nk >= 2, sum = w: L_{-n} X = :(d^{n-2}T) X:/(n-2)!"""
    def parts(n, mx):
        if n == 0: yield (); return
        for k in range(min(n, mx), 1, -1):
            for rest in parts(n - k, k): yield (k,) + rest
    out = []
    for pt in parts(w, w):
        X = {(): F(1)}
        for n in reversed(pt):
            dT = T
            for _ in range(n - 2): dT = K.dpoly(dT)
            X = K.scale(K.ope(dT, X).get(0, {}), F(1, factorial(n - 2))) if X != {(): F(1)} else K.scale(dT, F(1, factorial(n - 2)))
        out.append((pt, X))
    return out
def in_span(fields, G):
    e = K.Ech()
    for f_ in fields: e.insert(dict(f_))
    return not e.reduce(dict(G))
def field_conserved(G, alpha):
    """field level: Res(G, alpha) in the image of D_alpha (kwick's commutant test for one field)"""
    W = K.weight(next(iter(G))); qm = K.monomials(W - 2)
    e = K.Ech()
    for q in qm: e.insert(K.D_exp(q, alpha))
    return not e.reduce(K.residue_exp(G, alpha))
# ---------------- (i') the same cylinder rule by a memoised recursion per component (for large fields); must agree with eigen_raw
_EIG = {}
def _eig_letters(letters):
    """letters: sorted tuple of (c, q) of ONE component; returns the sympy polynomial in P (symbol PA) of the rule."""
    if not letters: return sp.Integer(1)
    key = (letters, K.KAPPA)
    if key in _EIG: return _EIG[key]
    kap = sp.Rational(K.KAPPA.numerator, K.KAPPA.denominator)
    first, rest = letters[0], letters[1:]
    tot = 0
    if first[1] == 0: tot += PA * _eig_letters(rest)
    for i, o in enumerate(rest):
        tot += kap * (-1)**o[1] * RD[first[1] + o[1]] * _eig_letters(rest[:i] + rest[i + 1:])
    tot = sp.expand(tot); _EIG[key] = tot; return tot
def eigen_fast(G):
    tot = 0
    for m, c in G.items():
        mv = tuple(l for l in m if l[0] == A); mz = tuple(l for l in m if l[0] == B)
        tot += sp.Rational(c.numerator, c.denominator) * _eig_letters(mv) * _eig_letters(mz).subs(PA, PB)
    return sp.expand(tot)
def solve_unique_Q(vectors, rhs):
    """x with sum_j x_j v_j = rhs when vectors are independent; None if inconsistent."""
    M, keys = build_mat(vectors + [rhs]); n = len(vectors)
    R, rk = M.rref()
    piv = [next(j for j in range(n + 1) if R[i, j] != 0) for i in range(rk)]
    if n in piv or rk != n: return None
    return [F(int(R[i, n].p), int(R[i, n].q)) for i in range(n)]
def exp_pole(G, alpha, r):
    """coefficient of (z-w)^{-r} in G(z) e^{alpha.Phi}(w) divided by the exponential (r = 1: kwick.residue_exp)."""
    from itertools import combinations
    out = {}
    for m, c in G.items():
        n = len(m)
        for k in range(0, n + 1):
            for S in combinations(range(n), k):
                fac = F(c); P = 0
                for i in S:
                    cc, q = m[i]; fac *= K.KAPPA * alpha[cc] * ((-1) ** q * factorial(q)); P += q + 1
                if P - r < 0: continue
                rem = [m[i] for i in range(n) if i not in S]
                for sh, cf in Wk._taylor(rem, P - r):
                    nm = tuple(sorted(sh)); out[nm] = out.get(nm, 0) + fac * cf
    return {k: v for k, v in out.items() if v}
