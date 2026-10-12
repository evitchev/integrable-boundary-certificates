"""VIR1 details (post-hoc, labelled): the symbolic I_5(t) per solution, the
limit classes at the two-dimensional kernels, and small facts used in NOTE.md.
Exit 0 iff all assertions hold."""
import sys, time
from fractions import Fraction as F
import sympy as sp
sys.dont_write_bytecode = True
import vir_lib as V

t = sp.Symbol('t')
t0_ = time.time()
C6 = V.Commutant(6)
B6 = C6.B
i6 = {m: i for i, m in enumerate(B6)}
f6, _, _ = V.quotient_functionals(6, V.hv_span(6))
fails = 0

# (b) largest loop power in the weight-4 x weight-6 residues
mx = 0
for p in V.basis(4):
    for q in B6:
        for mono, powd in V.ope(p, q, -1).items():
            mx = max(mx, max(powd))
print('largest loop power in the weight-4 x weight-6 residues: %d' % mx)
fails += mx > 2

Nt = (t**2 - 25) / (t**2 - 1)
st = -24 * t / (t**2 - 1)
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

def R(x):
    return sp.Rational(x.numerator, x.denominator)

sym = {}
for sol in (1, 2, 3):
    P = P4sym(sol)
    # check the symbolic P4 against the Fraction one
    Nq, sq = V.curve_point(F(7, 2))
    num = V.P4(sol, Nq, sq)
    assert all(sp.cancel(P[m].subs(t, sp.Rational(7, 2)) - R(num[m])) == 0 for m in num)
    cols = []
    for i in C6.comp:
        vec = [sp.Integer(0)] * len(C6.Bt)
        for pm, pc in P.items():
            for mono, powd in V.ope(pm, B6[i], -1).items():
                vec[C6.it[mono]] += pc * sum(R(c) * Nt**k for k, c in powd.items())
        cols.append([sp.cancel(sum(R(a) * b for a, b in zip(lam, vec) if a and b != 0)) for lam in C6.proj])
    M = sp.Matrix(len(C6.proj), len(C6.comp), lambda r, c: cols[c][r])
    from sympy.polys.matrices import DomainMatrix
    Mp = (M * (t**2 - 1)**4).applyfunc(sp.cancel)
    dM = DomainMatrix.from_Matrix(Mp).convert_to(sp.QQ[t]).to_field()
    nsd = dM.nullspace().to_Matrix()
    assert nsd.rows == 1, nsd.shape
    v = [sp.cancel(x) for x in nsd.row(0)]
    assert all(sp.cancel(e) == 0 for e in (M * sp.Matrix(v))), 'symbolic kernel vector does not solve the system'
    den = sp.lcm([sp.fraction(x)[1] for x in v])
    v = [sp.cancel(x * den) for x in v]
    g = sp.gcd_list(v)
    v = [sp.cancel(x / g) for x in v]
    deg = max(sp.Poly(x, t).degree() for x in v if x != 0)
    sym[sol] = v
    # non-Virasoro coordinates, symbolically
    full = [sp.Integer(0)] * len(B6)
    for c, i in zip(v, C6.comp):
        full[i] = c
    coords = [sp.expand(sum(R(a) * b for a, b in zip(lam, full) if a and b != 0)) for lam in f6]
    xtop = full[i6[((1, 1, 1, 1, 1, 1), ())]]
    print('Sol %d: symbolic I_5(t): primitive polynomial vector of degree %d; non-Virasoro coordinates %s; '
          '(dX)^6 coefficient %s   (%.0fs)' % (sol, deg, coords, sp.factor(xtop), time.time() - t0_), flush=True)
    fails += any(c != 0 for c in coords)
    sym[sol] = full

def limit_class(sol, t0):
    return [F(int(x.p), int(x.q)) for x in [sp.Rational(sp.cancel(c.subs(t, R(t0)))) for c in sym[sol]]]

print()
for name, t0, sols in (('jump point t = 5/3 (N = -25/2)', F(5, 3), (2, 3)), ('anchor t = 3 (N = -2)', F(3), (1, 2, 3))):
    N, s = V.curve_point(t0)
    lims = {sol: limit_class(sol, t0) for sol in sols}
    for sol in sols:
        K = C6.kernel(V.P4(sol, N, s), N)
        inK = V.rank(K + [lims[sol]]) == len(K)
        nz = any(lims[sol])
        print('%s: Sol %d limit class nonzero %s, lies in the kernel of ITS OWN I_3 there (dim %d): %s; non-Virasoro coords %s'
              % (name, sol, nz, len(K), inK, [str(c) for c in V.apply(f6, lims[sol])]))
        fails += (not inK) or (not nz)
    same_P4 = [(a, b, V.P4(a, N, s) == V.P4(b, N, s)) for a in sols for b in sols if a < b]
    print('   same I_3 density at the point: %s' % same_P4)
    for a in sols:
        for b in sols:
            if a < b and V.P4(a, N, s) == V.P4(b, N, s):
                K = C6.kernel(V.P4(a, N, s), N)
                r = V.rank([lims[a], lims[b]])
                print('   Sols %d,%d share I_3: their two limit classes have rank %d; kernel dim %d -> they %s the kernel'
                      % (a, b, r, len(K), 'SPAN' if r == len(K) else 'do NOT span'))
# the anchor kernel of Sol 3: pure-Y and pure-X directions
N, s = V.curve_point(3)
K = C6.kernel(V.P4(3, N, s), N)
ypure = [i for i, m in enumerate(B6) if m[0] == ()]
xpure = [i for i, m in enumerate(B6) if m[1] == ()]
mixed = [i for i, m in enumerate(B6) if m[0] != () and m[1] != ()]
rk_mixed = V.rank([[v[i] for i in mixed] for v in K]) if any(v[i] for v in K for i in mixed) else 0
print('anchor, Sol 3: kernel dim %d; rank of the mixed (X and Y) components over the kernel: %d '
      '(0 = the kernel is {pure-X class} + {pure-Y class})' % (len(K), rk_mixed))
print('I_3 at the anchor (Sol 3):', {k: str(v) for k, v in V.P4(3, N, s).items()})
print('DETAILS', 'FAILED' if fails else 'OK')
sys.exit(1 if fails else 0)
