"""VIR2: symbolic I_5(t) per solution (primitive polynomial vectors over V.basis(6)),
same method as VIR1 v4_details.py; saved to i5_sym.pkl.  Exit 0 iff verified."""
import sys, time, pickle
from fractions import Fraction as F
import sympy as sp
from sympy.polys.matrices import DomainMatrix
sys.dont_write_bytecode = True
import vir_lib as V
t = sp.Symbol('t')
C6 = V.Commutant(6)
B6 = C6.B
Nt = (t**2 - 25) / (t**2 - 1)
st = -24 * t / (t**2 - 1)
def R(x):
    return sp.Rational(x.numerator, x.denominator)
def P4sym(sol):
    N, s = Nt, st
    base = {((), ((1, 1), (1, 1))): sp.Integer(1), ((), ((2, 2),)): sp.Integer(-2)}
    if sol == 1:
        ex = {((2, 2), ()): (-9 - 3*N - s)/6, ((1, 1), ((1, 1),)): sp.Integer(6), ((1, 1, 1, 1), ()): (2 + N + s)/3}
    elif sol == 2:
        ex = {((2, 2), ()): -(409 + 46*N + N*N + (N - 5)*s)/192, ((1, 1), ((1, 1),)): (11 + N + s)/8,
              ((1, 1, 1, 1), ()): (117 + 2*N + N*N + (15 + N)*s)/192}
    else:
        ex = {((2, 2), ()): (13 - 27*N + 2*N*N + (5 - 2*N)*s)/6, ((1, 1), ((1, 1),)): 7 - N + s, ((1, 1, 1, 1), ()): sp.Integer(1)}
    base.update(ex)
    return {m: sp.cancel(c) for m, c in base.items()}
out = {}
t0 = time.time()
ok = True
for sol in (1, 2, 3):
    P = P4sym(sol)
    cols = []
    for i in C6.comp:
        vec = [sp.Integer(0)] * len(C6.Bt)
        for pm, pc in P.items():
            for mono, powd in V.ope(pm, B6[i], -1).items():
                vec[C6.it[mono]] += pc * sum(R(c) * Nt**k for k, c in powd.items())
        cols.append([sp.cancel(sum(R(a_) * b_ for a_, b_ in zip(lam, vec) if a_ and b_ != 0)) for lam in C6.proj])
    M = sp.Matrix(len(C6.proj), len(C6.comp), lambda r, c: cols[c][r])
    Mp = (M * (t**2 - 1)**4).applyfunc(sp.cancel)
    nsd = DomainMatrix.from_Matrix(Mp).convert_to(sp.QQ[t]).to_field().nullspace().to_Matrix()
    assert nsd.rows == 1
    v = [sp.cancel(x) for x in nsd.row(0)]
    assert all(sp.cancel(e) == 0 for e in (M * sp.Matrix(v)))
    den = sp.lcm([sp.fraction(x)[1] for x in v])
    v = [sp.cancel(x * den) for x in v]
    gg = sp.gcd_list(v)
    v = [sp.expand(sp.cancel(x / gg)) for x in v]
    full = [sp.Integer(0)] * len(B6)
    for c, i in zip(v, C6.comp):
        full[i] = c
    # cross-check against the numeric kernel at two fibres
    for tv in (F(7, 2), F(-9, 4)):
        N, s = V.curve_point(tv)
        K = C6.kernel(V.P4(sol, N, s), N)
        num = [F(int(sp.Rational(c.subs(t, R(tv))).p), int(sp.Rational(c.subs(t, R(tv))).q)) for c in full]
        if not (len(K) == 1 and V.rank([K[0], num]) == 1 and any(num)):
            ok = False
    out[sol] = full
    print('Sol %d: degree %d, verified vs numeric kernel: %s  (%.0fs)' % (sol, max(sp.Poly(x, t).degree() for x in full if x != 0), ok, time.time() - t0), flush=True)
pickle.dump({k: [sp.srepr(x) for x in v] for k, v in out.items()}, open('i5_sym.pkl', 'wb'))
sys.exit(0 if ok else 1)
