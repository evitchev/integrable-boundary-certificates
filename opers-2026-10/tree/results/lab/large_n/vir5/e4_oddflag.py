"""VIR5a: the odd-order flag at k = -6 (t = 4).  The flag in the prediction file came from the
Beta-ratio formula, whose Pochhammer (i-1)/n + j hits zero for odd i when n = -2.  Direct test,
independent of that formula: is the odd-order WKB term an exact derivative, T_i = theta(Q) with
Q in the ring?  (Then its integral vanishes and there is no even-spin charge.)
Also run at k = 10/3 as a control where the flag was True.  Exit 0 iff exact at every odd order."""
import sys
from fractions import Fraction as F
import sympy as sp
sys.dont_write_bytecode = True
import mellin_wkb as MW

def exact(eng, i):
    """Solve theta(Q) = T0 * W_i in the ring; return (is_exact, note)."""
    n = eng.n
    Ti = MW.el_mul({(1, 1 / n): MW.ONE}, eng.W[i])
    avals = {a for (a, b) in Ti}
    assert len(avals) == 1
    a = avals.pop()
    bs = sorted({b for (_, b) in Ti}, reverse=True)
    # all b's differ by integers
    b = bs[0]; q_prev = MW.ZERO; note = ''
    while b >= bs[-1] - 1:
        c = Ti.get((a, b), MW.ZERO)
        # coefficient at p^b:  q_b (a + hM b) - hM (b+1) q_{b+1} = c_b
        rhs = c + q_prev * MW.P_(MW.R_(eng.hM * (b + 1)))
        fac = MW.R_(F(a) + eng.hM * b)
        if b < bs[-1]:
            return rhs.is_zero, note
        if fac == 0:
            if not rhs.is_zero:
                return False, 'resonant obstruction at p^%s' % b
            q = MW.ZERO; note = 'resonance at p^%s (rhs zero)' % b
        else:
            q = rhs * MW.P_(1 / fac)
        q_prev = q
        b -= 1
    return None, 'not reached'

fails = 0
for k in (sp.Rational(-6), sp.Rational(10, 3), sp.Rational(1, 2)):
    n = k + 4
    eng = MW.MellinWKB(n, 1 / k, MW.family_symbol(k, 6))
    eng.solve(9)
    out = {i: exact(eng, i) for i in (3, 5, 7, 9)}
    print('k = %s (n = %s): odd-order terms exact derivatives: %s' % (k, n, out))
    fails += sum(1 for v in out.values() if v[0] is not True)
print('ODD ORDERS', 'all exact derivatives: no even-spin charge at any of these fibres' if not fails else 'NOT all exact (%d)' % fails)
sys.exit(1 if fails else 0)
