"""Kit verifier (v1.2.5, new; pure Python/sympy, no Singular or msolve) of EXC4 part B's archived solver output (Sol 3 level two, four apparent
singularities).  Run with cwd = results/lab/large_n/exc4.  Usage: gb_verify.py <TAG>   (TAG = t2_P1-2_Q1-3 | t9_4_P1-2_Q1-3 | t2_P1-3_Q1-5 | t9_4_P2-5_Q1-4)
Inputs (all archived): the polynomial system h_fourq_<TAG>_four.sing (written by h_four_q.py; ideal I) and its msolve input h_fourq_<TAG>_four.ms;
msolve's grevlex Groebner basis G as passed to Singular in gb_elim_<TAG>.sing; Singular's fglm eliminants gb_elim_<TAG>.out; the data block
d2_sol3_<t>_P<P>_Q<Q>.json (record's tool); PREDICTION_EXC4.json (hold-out points).
Checks, all in exact rational arithmetic:
 (1) the msolve input .ms is the same polynomial set as the Singular input ideal I (up to sign);
 (2) G is reduced and its leading monomials leave exactly 40 standard monomials (vdim 40); the multiplication matrices of all 12 variables
     (normal forms computed from G) COMMUTE, which certifies G as a Groebner basis of the ideal it generates (so they are the true multiplication maps);
 (3) every input polynomial has normal form 0 modulo G (so every solution of G solves the system: V(G) is contained in V(I));
 (4) M_ii = I_1 * identity with I_1 = vacuum + 4 (level two);  quintic(M_ee) = 0 with the quintic irreducible, so charpoly(M_ee) = quintic^8 :
     the 40 ordered solutions carry exactly the five I_3 values, eight orderings each;  the characteristic polynomial of the linear form
     z + 3 z3 + 5 r12 + 7 u12 is squarefree of degree 40 (40 distinct ordered solutions = 5 unordered configurations x 8 relabellings);
 (5) the quintic equals Singular's archived fglm eliminant, the registered prediction (hold-out points) and the data block's quintic factor.
Not checked here: that G generates ALL of I (no solution lost; that is msolve's computation) -- (1)-(4) show every solution of G is a solution of I."""
import json, re, sys, time
from fractions import Fraction
import sympy as sp
from sympy.polys.rings import ring
from sympy.polys.domains import QQ
from sympy.polys.orderings import grevlex
from sympy.polys.matrices import DomainMatrix
TAG = sys.argv[1]; mode = sys.argv[2] if len(sys.argv) > 2 else 'real'
T0 = time.time(); ok = True
def say(*a): print(*a, flush=True)
def block(txt, start):
    i = txt.index(start) + len(start); j = txt.index(';', i); return [p.strip() for p in re.split(r',\s*\n', txt[i:j]) if p.strip()]
sing_in = open('h_fourq_%s_four.sing' % TAG).read()
VARS = re.search(r'ring R = 0,\(([^)]*)\),dp;', sing_in).group(1).split(',')
gbs = open('gb_elim_%s.sing' % TAG).read()
assert re.search(r'ring R = 0,\(([^)]*)\),dp;', gbs).group(1).split(',') == VARS
R, *X = ring(','.join(VARS), QQ, grevlex); V = dict(zip(VARS, X)); nv = len(VARS)
def parse(s):
    return R(sp.sympify(s.replace('^', '**'), locals={v: sp.Symbol(v) for v in VARS}))
I = [parse(p) for p in block(sing_in, 'ideal I =')]
ms = open('h_fourq_%s_four.ms' % TAG).read().split('\n', 2)
assert ms[0].replace(' ', '').split(',') == VARS and ms[1].strip() == '0'
MS = [parse(p) for p in re.split(r',\s*\n', ms[2].strip().rstrip(',')) if p.strip()]
G = [parse(p) for p in block(gbs, 'ideal G =')]
mon = lambda f: f.monic()
same_set = sorted(str(mon(f)) for f in I) == sorted(str(mon(f)) for f in MS)
say('%s: %d input polynomials (Singular .sing), %d (msolve .ms), same set: %s; Groebner basis G: %d elements (%.0fs)' % (TAG, len(I), len(MS), same_set, len(G), time.time() - T0))
ok &= same_set
if mode == 'tamper':    # negative control: perturb one input polynomial; it must no longer reduce to 0
    I[0] = I[0] + R(1)
LM = {g.LM: g for g in G}
def divisible(m): return any(all(a >= b for a, b in zip(m, l)) for l in LM)
std = []; frontier = [tuple([0] * nv)]; seen = set(frontier)
while frontier:
    m = frontier.pop()
    if divisible(m): continue
    std.append(m)
    for i in range(nv):
        m2 = tuple(m[j] + (1 if j == i else 0) for j in range(nv))
        if m2 not in seen and sum(m2) < 60: seen.add(m2); frontier.append(m2)
std.sort(); idx = {m: i for i, m in enumerate(std)}; N = len(std)
reduced = all(not divisible(m) for g in G for m, c in g.terms() if m != g.LM)
say('   standard monomials (vdim): %d; G is a reduced basis (no tail term divisible by a leading monomial): %s (%.0fs)' % (N, reduced, time.time() - T0))
ok &= N == 40 and reduced
# normal forms of monomials as vectors on the standard monomials (recursive, memoised; valid for a reduced Groebner basis, certified below)
NF = {}
def nf(u):
    if u in NF: return NF[u]
    if u in idx: v = {idx[u]: QQ(1)}
    elif u in LM:
        g = LM[u]; lc = g.LC; v = {}
        for m, c in g.terms():
            if m != u: v[idx[m]] = -c / lc
    else:
        l = next(l for l in LM if all(a >= b for a, b in zip(u, l)))
        k = next(j for j in range(nv) if u[j] > l[j])
        w = nf(tuple(u[j] - (j == k) for j in range(nv))); v = {}
        for si, c in w.items():
            for ti, d in nf(tuple(std[si][j] + (j == k) for j in range(nv))).items():
                v[ti] = v.get(ti, QQ(0)) + c * d
        v = {i: c for i, c in v.items() if c != 0}
    NF[u] = v; return v
def nf_poly(f):
    v = {}
    for m, c in f.terms():
        for i, d in nf(m).items(): v[i] = v.get(i, QQ(0)) + c * d
    return {i: c for i, c in v.items() if c != 0}
def mult_matrix(j):
    rows = [[QQ(0)] * N for _ in range(N)]
    for si, m in enumerate(std):
        for ti, c in nf(tuple(m[q] + (q == j) for q in range(nv))).items(): rows[ti][si] = c
    return DomainMatrix(rows, (N, N), QQ)
M = {v: mult_matrix(j) for j, v in enumerate(VARS)}
say('   multiplication matrices of the %d variables built (%.0fs)' % (nv, time.time() - T0))
comm = all((M[a] * M[b]).to_Matrix() == (M[b] * M[a]).to_Matrix() for i, a in enumerate(VARS) for b in VARS[i + 1:])
say('   all %d pairs of multiplication matrices commute (G is a Groebner basis): %s (%.0fs)' % (nv * (nv - 1) // 2, comm, time.time() - T0)); ok &= comm
nz = sum(1 for f in I if nf_poly(f))
say('   input polynomials with normal form 0 modulo G: %d of %d (every solution of G solves the system) (%.0fs)' % (len(I) - nz, len(I), time.time() - T0)); ok &= nz == 0
e = sp.Symbol('e')
outt = open('gb_elim_%s.out' % TAG).read()
pe = sp.Poly(sp.sympify(re.search(r'eliminant ee:\n(.*)', outt).group(1).replace('^', '**'), locals={'ee': e}), e).monic()
pi = sp.Poly(sp.sympify(re.search(r'eliminant ii:\n(.*)', outt).group(1).replace('^', '**'), locals={'ii': e}), e).monic()
assert pi.degree() == 1
I1 = -pi.all_coeffs()[1]
dense = lambda f: DomainMatrix([[f(i, j) for j in range(N)] for i in range(N)], (N, N), QQ)
Id = dense(lambda i, j: QQ(1) if i == j else QQ(0)); Z0 = dense(lambda i, j: QQ(0))
ii_scalar = M['ii'].to_Matrix() == (Id * QQ(I1.p, I1.q)).to_Matrix()
quint = pe
if mode == 'tamper_quintic': quint = sp.Poly(pe.as_expr() + sp.Rational(1, 1000), e)
acc = Z0
for c in quint.all_coeffs(): acc = acc * M['ee'] + Id * QQ(sp.Rational(c).p, sp.Rational(c).q)
annih = acc.to_Matrix() == Z0.to_Matrix()
irred = len(sp.factor_list(quint.as_expr())[1]) == 1 and quint.degree() == 5
say('   M_ii = %s * identity: %s;  quintic(M_ee) = 0: %s;  quintic irreducible over Q: %s  =>  charpoly(M_ee) = quintic^8 (%.0fs)' % (I1, ii_scalar, annih, irred, time.time() - T0))
ok &= ii_scalar and annih and irred
Ml = M['z'] + M['z3'] * QQ(3) + M['r12'] * QQ(5) + M['u12'] * QQ(7)   # a linear form separating the ordered solutions
cz = Ml.charpoly(); cz = sp.Poly([sp.Rational(int(c.numerator), int(c.denominator)) for c in cz], e)
sqf = sp.gcd(cz, cz.diff(e)).degree() == 0
say('   charpoly of M_(z + 3 z3 + 5 r12 + 7 u12): degree %d, squarefree: %s (40 distinct ordered solutions) (%.0fs)' % (cz.degree(), sqf, time.time() - T0)); ok &= sqf and cz.degree() == 40
m = re.match(r't(\d+(?:_\d+)?)_P(\d+)-(\d+)_Q(\d+)-(\d+)$', TAG)
tstr, Pn, Pd, Qn, Qd = m.groups()
d2 = json.load(open('d2_sol3_t%s_P%s_%s_Q%s_%s.json' % (tstr, Pn, Pd, Qn, Qd)))
facs = [sp.Poly(sp.sympify(f, locals={'e': e}), e) for f in d2['factors']]
dq = [f for f in facs if f.degree() == 5][0].monic(); lin = [f for f in facs if f.degree() == 1]
eq_data = dq == quint
say('   quintic == data block quintic factor (d2_sol3_*.json, factors %s): %s;  sixth state linear (level shift): %s' % ([f.degree() for f in facs], eq_data, len(lin) == 1))
ok &= eq_data and len(lin) == 1
P = json.load(open('PREDICTION_EXC4.json'))
reg = [p for p in P['points'] if p['t'] == tstr.replace('_', '/') and p['P'] == '%s/%s' % (Pn, Pd) and p['Q'] == '%s/%s' % (Qn, Qd)]
if reg:
    rq = sp.Poly(sp.sympify(reg[0]['I3_quintic_monic'], locals={'e': e}), e); ri = sp.Rational(reg[0]['I1_monic'])
    eq_reg = rq == quint and ri == I1 and reg[0]['vdim_ordered'] == N
    say('   REGISTERED hold-out point: registered quintic, I_1 and vdim equal the verified ones: %s' % eq_reg); ok &= eq_reg
else:
    say('   not a registered point (post hoc: t = 2, (1/2,1/3))')
say('GB VERIFY %s [%s]: %s (%.0fs)' % (TAG, mode, 'PASS' if ok else 'FAIL', time.time() - T0))
sys.exit(0 if ok else 1)
