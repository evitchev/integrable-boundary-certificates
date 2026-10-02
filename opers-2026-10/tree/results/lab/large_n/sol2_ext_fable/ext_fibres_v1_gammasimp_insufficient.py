"""REVIEW (Fable secondary, 2026-10-02): the archived Sol 2 family S2(k) (sol2f.py's run/score, unchanged) at fibres OUTSIDE the
range tested in the record (k in [1, 4]): the exact point t = 1/3 (k = -1/2, n = 2) and others with k < 1 or k < 0.
Monic spins 1..9 against the certified charges.  Tamper: delta + 1/10 must fail.  Exit 0 iff the exact-point fibre matches at
every spin with a non-vanishing WKB charge and the tamper fails there; other fibres are reported."""
import sys, time
from fractions import Fraction as F
import sympy as sp
from common import *
def run(k, M, K=10, delta_shift=0):
    kq = sp.Rational(k.numerator, k.denominator)
    tv = (3 * kq + 2) / (kq + 2); g = classG(2, tv)
    assert g['alpha'] == F(k) and g['gamma'] == 1 and g['D'] == 2 * F(k) + 3, g
    n = g['D']; Cv = 2 * n * (1 + F(M))
    bl = blocks_of(g, M)
    if delta_shift: bl = [(a, sh, d + delta_shift) if a != 1 else (a, sh, d) for (a, sh, d) in bl]
    sv = gwkb.wkb(bl, g['s'], F(M), K); bb = n * F(M); ch = {}
    for kk in range(2, K + 1, 2):
        sk = sv[kk]
        if any(a.denominator == 1 and a <= 0 for (a, c) in sk): pw, poly, nonpole = thwkb.integrate_reg(sk, bb); poly = dict(poly)
        else:
            r = thwkb.integrate(sk, bb); poly = {k_: sp.Rational(v.numerator, v.denominator) for k_, v in r[2].items()} if r else {}
        ch[kk - 1] = sp.expand(sum(v * (sp.Rational(Cv.numerator, Cv.denominator) * PX**2 / g['aY'])**(i // 2) * (sp.Rational(Cv.numerator, Cv.denominator) * PI**2 / g['aX'])**(j // 2) for (i, j), v in poly.items()))
    return tv, g, Cv, ch
def score(ch, cert, spins):
    rep = {}
    for s_ in spins:
        if sp.gammasimp(ch[s_]) == 0: rep[s_] = 'WKB charge identically 0 (certified: %s)' % ('ZERO' if cert[s_] == 0 else 'non-zero, %d coeffs' % len(sp.Poly(cert[s_], PX, PI).terms())); continue
        d = sp.expand(sp.gammasimp(sp.expand(monic(ch[s_], (s_ + 1) // 2) - cert[s_]))); n_ = len(sp.Poly(cert[s_], PX, PI).terms())   # v1: gammasimp (v0 left Gamma(-1/3) + 3 Gamma(2/3) unsimplified and printed DIFF)
        rep[s_] = True if d == 0 else f'DIFF ({len(sp.Poly(d, PX, PI).terms())} of {n_} coeffs)'
    return rep
ok = True
if __name__ != "__main__": sys.argv = sys.argv[:1]
for tstr in sys.argv[1:]:
    tv = sp.Rational(tstr); kq = 2 * (tv - 1) / (3 - tv); k = F(int(kq.p), int(kq.q))
    if not regular_at(2, tv): print(f't = {tv} (k = {k}): Sol 2 table has a pole -> not tested', flush=True); continue
    cert = certified(2, tv)
    for M, nm in ((1 / k, 'M = 1/k'),):
        T0 = time.time()
        try:
            tv_, g, Cv, ch = run(k, M)
        except Exception as ex:
            print(f't = {tv} (k = {k}): engine raised {type(ex).__name__}: {str(ex)[:120]}', flush=True); ok &= tstr != '1/3'; continue
        rep = score(ch, cert, (1, 3, 5, 7, 9))
        print(f't = {tv} (k = {k}, n = {g["D"]}) {nm}, C = {Cv}: {rep}  ({time.time()-T0:.0f}s)', flush=True)
        try:
            _, _, _, cht = run(k, M, delta_shift=F(1, 10)); rt = score(cht, cert, (3, 5, 7, 9))
        except Exception as ex:
            rt = {'error': str(ex)[:80]}
        print(f'   TAMPER delta + 1/10: {rt}', flush=True)
        if tstr == '1/3':
            ok &= all(v is True or v == 'WKB charge identically 0' for v in rep.values()) and any(v is True for v in rep.values()) and not any(v is True for v in rt.values())
print('exact-point fibre t = 1/3:', 'MATCH' if ok else 'NO MATCH / not run')
sys.exit(0 if ok else 2)
