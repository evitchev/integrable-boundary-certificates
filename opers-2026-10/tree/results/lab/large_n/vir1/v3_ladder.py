"""VIR1 lead addition 1 (extended): the August ladder I_7, I_9, I_11 with
QUANTUM products, on the repo's promoted Q(t)-symbolic charges
(results/qt_promotion.json, read-only; A(N) + s B(N) on the even mixed basis).
Post-seal extension; spins >= 7 were outside the seal.

For each solution and sigma in (7, 9, 11): the coordinates of the promoted
density in Q_w = V_w / (HV_w + dV_{w-1}), w = sigma + 1,
  (i)  exactly at 24 curve points (Fractions);
  (ii) SYMBOLICALLY in N: sum_i lam_i A_i(N) == 0 and sum_i lam_i B_i(N) == 0.
Controls: classical ladder entries of the bcd note at t = 2; a plant; vacuity;
sigma 7 promoted density == my own kernel class at t = 2.
Exit 0 all in HV + d and controls fire; 2 some charge outside; 1 control failure.
"""
import json, sys, time, random
import os as _os_kit; _KITROOT = _os_kit.path.join(_os_kit.path.dirname(_os_kit.path.abspath(__file__)), '..', '..', '..', '..') + '/'  # KIT PATCH: repo root from __file__ (was an absolute path)
from fractions import Fraction as F
import sympy as sp
sys.dont_write_bytecode = True
sys.path.insert(0, _KITROOT + 'code')
import vir_lib as V
from mixed_engine import gen_mixed_basis

QT = json.load(open(_KITROOT + 'results/qt_promotion.json'))
Ns = sp.Symbol('N')
t0 = time.time()
LOOP0 = F(7, 3)
bad = ctlbad = 0
TS = [F(2), F(-2), F(7, 2), F(-7, 2), F(4), F(-4), F(11, 5), F(-11, 5), F(-3, 7), F(3, 7), F(9, 2), F(-9, 2),
      F(-13, 4), F(13, 4), F(1, 2), F(-1, 2), F(6), F(-6), F(8, 3), F(-8, 3), F(1, 5), F(-1, 5), F(10), F(-10)]
MAXSPIN = int(sys.argv[1]) if len(sys.argv) > 1 else 11

def frac(x):
    x = sp.nsimplify(x) if not x.is_Rational else x
    return F(int(x.p), int(x.q))

for sigma in (7, 9, 11):
    if sigma > MAXSPIN:
        break
    w = sigma + 1
    B = V.basis(w)
    assert B == gen_mixed_basis(w), 'basis order differs from the repo'
    idx = {m: i for i, m in enumerate(B)}
    hv = V.hv_span(w, LOOP0)
    funcs, _, _ = V.quotient_functionals(w, hv)
    print('\n== sigma %d (weight %d): dim V = %d, HV span %d, dim Q_w = %d  (%.0fs)'
          % (sigma, w, len(B), len(hv), len(funcs), time.time() - t0), flush=True)
    # vacuity
    rng = random.Random(100 + sigma)
    g = [F(rng.randint(-9, 9), rng.randint(1, 9)) for _ in B]
    nz = sum(1 for c in V.apply(funcs, g) if c)
    print('   vacuity: generic vector has %d of %d coordinates nonzero' % (nz, len(funcs)))
    if nz != len(funcs):
        ctlbad += 1
    drows = [V.to_row(V.d_mono(m), idx) for m in V.basis(w - 1)]
    def member(vec, span):
        rows = [V.to_row(v, idx) for v in span] + drows
        r0 = V.rank(rows)
        return V.rank(rows + [vec]) == r0
    for sol in (1, 2, 3):
        data = QT['cyl%d_sigma%d' % (sol, sigma)]
        assert data['symbolic_verified'] is True
        AB = {int(i): (sp.sympify(a, locals={'N': Ns}), sp.sympify(b, locals={'N': Ns}))
              for i, (a, b) in data['coeffs_A_B_of_N'].items()}
        # (ii) symbolic
        symzero = True
        for lam in funcs:
            SA = sp.cancel(sum(sp.Rational(c.numerator, c.denominator) * AB[i][0]
                               for i, c in enumerate(lam) if c and i in AB))
            SB = sp.cancel(sum(sp.Rational(c.numerator, c.denominator) * AB[i][1]
                               for i, c in enumerate(lam) if c and i in AB))
            if SA != 0 or SB != 0:
                symzero = False
        # (i) pointwise
        ptbad = 0
        skipped = 0
        vec2 = None
        for t in TS:
            N, s = V.curve_point(t)
            Nq, sq = sp.Rational(N.numerator, N.denominator), sp.Rational(s.numerator, s.denominator)
            vec = [F(0)] * len(B)
            pole = False
            for i, (a, b) in AB.items():
                val = a.subs(Ns, Nq) + sq * b.subs(Ns, Nq)
                if not val.is_Rational:
                    pole = True
                    break
                vec[i] = F(int(val.p), int(val.q))
            if pole:
                skipped += 1
                continue
            if any(V.apply(funcs, vec)):
                ptbad += 1
            if t == 2:
                vec2 = vec
        print('   Sol %d I_%d: symbolic in N: all %d coordinates == 0: %s; pointwise outside at %d of %d fibres (%d skipped)  (%.0fs)'
              % (sol, sigma, len(funcs), symzero, ptbad, len(TS) - skipped, skipped, time.time() - t0), flush=True)
        if not symzero or ptbad:
            bad += 1
        # controls at t = 2
        N, s = V.curve_point(2)
        pl = list(vec2); pl[idx[((), ((w // 2, w // 2),))]] += F(1, 1000)
        if not any(V.apply(funcs, pl)):
            ctlbad += 1; print('   CONTROL-FAIL plant not flagged')
        if sigma == 7:
            C8 = V.Commutant(8)
            K = C8.kernel(V.P4(sol, N, s), N)
            same = len(K) == 1 and V.rank(drows + [K[0], vec2]) == V.rank(drows) + 1
            print('      promoted I_7 at t=2 == my own kernel class mod d: %s' % same)
            if not same:
                ctlbad += 1
        if sigma in (7, 9):
            cl_lo = member(vec2, V.classical_span(w, w - 2))
            print('      classical at t=2: I_%d in A<=%d: %s (record: no)' % (sigma, w - 2, cl_lo), flush=True)
            if cl_lo:
                ctlbad += 1
print()
print('RESULT:', 'CONTROL FAILURE' if ctlbad else ('%d charges OUTSIDE <X,T_Y>' % bad if bad else
      'every promoted charge tested lies in the vertex subalgebra <X, T_Y> + Im d, identically in N'))
sys.exit(1 if ctlbad else (2 if bad else 0))
