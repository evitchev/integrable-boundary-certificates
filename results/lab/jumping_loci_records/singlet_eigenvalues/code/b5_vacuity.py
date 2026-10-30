"""B5 (iii): rank of the cylinder vacuum-eigenvalue map on one-boson fields of weight w (the records' VIR3-R recipe: pairings with
C(a,b) = (-1)^b R^(a+b)(0), R(eps) = 1/(4 sinh^2(eps/2)) - 1/eps^2; unpaired letters must be dV, each giving P). Exact Fractions; kappa = 1 (rank is kappa-independent)."""
import sys, hashlib, os
from fractions import Fraction as F
from math import factorial
import sympy as sp
print('script sha256', hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest())
eps = sp.Symbol('eps')
Rser = sp.series(1/(4*sp.sinh(eps/2)**2) - 1/eps**2, eps, 0, 24).removeO()
Rd = []
for k in range(20):
    v = sp.Rational(sp.diff(Rser, eps, k).subs(eps, 0)); Rd.append(F(int(v.p), int(v.q)))
print('R^(k)(0), k=0..8:', Rd[:9])
def C(a, b): return (-1)**b * Rd[a + b]
def partitions(n, mx=None):
    if mx is None: mx = n
    if n == 0: yield (); return
    for k in range(min(n, mx), 0, -1):
        for rest in partitions(n - k, k): yield (k,) + rest
def pairings(L):
    if not L: yield [], []; return
    first, rest = L[0], L[1:]
    for pr, un in pairings(rest): yield pr, [first] + un
    for i, o in enumerate(rest):
        for pr, un in pairings(rest[:i] + rest[i+1:]): yield [(first, o)] + pr, un
def eigen(mono):  # mono: tuple of derivative orders k (letter q = k-1); returns dict degree -> coefficient
    qs = [k - 1 for k in mono]; out = {}
    for pr, un in pairings(qs):
        if any(qq >= 1 for qq in un): continue
        val = F(1)
        for (a, b) in pr: val *= C(a, b)
        out[len(un)] = out.get(len(un), 0) + val
    return {d: v for d, v in out.items() if v}
def rank(rows):
    M = [list(r) for r in rows]; rk = 0; ncol = len(M[0]) if M else 0
    for c in range(ncol):
        piv = next((i for i in range(rk, len(M)) if M[i][c] != 0), None)
        if piv is None: continue
        M[rk], M[piv] = M[piv], M[rk]
        for i in range(len(M)):
            if i != rk and M[i][c] != 0:
                f = M[i][c] / M[rk][c]; M[i] = [x - f*y for x, y in zip(M[i], M[rk])]
        rk += 1
    return rk
def p(n): return sum(1 for _ in partitions(n))
for w in (3, 5, 7, 9):
    monos = list(partitions(w)); degs = list(range(w, 0, -2))
    rows = [[eigen(m).get(d, F(0)) for d in degs] for m in monos]
    rk = rank(rows); shapes = len(degs); modd = p(w) - p(w - 1)
    print(f'weight {w}: {len(monos)} monomials; classes mod d = p(w)-p(w-1) = {modd}; odd-polynomial shape space dim = {shapes}; rank of the eigenvalue map = {rk}')
    if w <= 7:
        for m, r in zip(monos, rows): print(f'   {m}: ' + ' + '.join(f'({c})P^{d}' for c, d in zip(r, degs) if c) if any(r) else f'   {m}: 0')
    verdict = ('AUTOMATIC: every weight-%d one-boson eigenvalue is a multiple of one polynomial' % w) if rk == 1 else (('NON-AUTOMATIC: rank %d < shape dim %d, a %d-parameter family up to scale' % (rk, shapes, rk - 1)) if rk < shapes else ('VACUOUS AS A DISCRIMINATOR: rank %d = full odd shape space (every odd degree-%d polynomial is some field\'s eigenvalue)' % (rk, w)))
    print('   ->', verdict)
# WUW1 cross-check: W_U = (4/3)(dU)^3 + 2 dU d^2U + (1/3) d^3U -> eigen = (4/3)(P^3 - P/4)
e = {}
for coef, mono in ((F(4,3), (1,1,1)), (F(2), (2,1)), (F(1,3), (3,))):
    for d, v in eigen(mono).items(): e[d] = e.get(d, 0) + coef*v
print('WUW1 cross-check eigen(W_U) =', {d: v for d, v in e.items() if v}, '-> expected {3: 4/3, 1: -1/3}:', 'PASS' if e == {3: F(4,3), 1: F(-1,3)} else 'FAIL')
sys.exit(0 if e == {3: F(4,3), 1: F(-1,3)} else 1)
