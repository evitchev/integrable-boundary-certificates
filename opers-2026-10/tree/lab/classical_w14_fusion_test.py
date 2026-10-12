"""Test of the fusion-product prediction at q = 8 (weight 14, t = 9/7, sheet 3): the extra invariant from the fast-path
witness record, its vacuum eigenvalue by the generating-function recipe, against -prod_{m=-3}^{3} Lambda(P + i m/2),
Lambda = P^2 - Q^2 + 49/16."""
import sys, os, json, time
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); os.chdir(os.path.join(ROOT, 'code')); sys.path.insert(0, '.'); sys.path.insert(0, os.path.join(ROOT, 'lab'))
from fractions import Fraction as F
import sympy as sp
from math import comb
from mixed_engine import gen_mixed_basis
import vev_genfun as G
W = 14; q = 8
rec = json.load(open(os.path.join(ROOT, 'results', 'lab', 'fast_kernel2', 'w14_t97_sheet3.json')))
BASIS = gen_mixed_basis(W); print("record:", rec['status'], "genuine_dim", rec['genuine_dim'], "; basis size", len(BASIS))
def excess(mono):
    xs, yp = mono; return sum(k - 1 for k in xs) + sum(k + l - 2 for k, l in yp)
def top_monomial(i): return ((1,) * (W - 2 * i), tuple(((1, 1),) * i))
ker = []
for w in rec['witnesses']:
    v = {BASIS[int(i)]: F(c) for i, c in w['charge_v'].items()}
    ker.append({m: c for m, c in v.items() if c})
k = W // 2; tops = [{m: c for m, c in v.items() if excess(m) == 0} for v in ker]
lam = sp.symbols(f'l0:{len(ker)}')
target = {top_monomial(i): F(comb(k, i) * (-1)**i) for i in range(k + 1)}
eqs = [sum(l * sp.Rational(tp.get(top_monomial(i), F(0)).numerator, tp.get(top_monomial(i), F(0)).denominator) for l, tp in zip(lam, tops)) - sp.Rational(target[top_monomial(i)].numerator, target[top_monomial(i)].denominator) for i in range(k + 1)]
so = sp.solve(eqs, lam, dict=True); print("combination with top part (x^2 - r^2)^7:", "found" if so else "NONE")
coef = [F(int(sp.Rational(so[0].get(l, 0)).p), int(sp.Rational(so[0].get(l, 0)).q)) for l in lam]
E = {}
for c_, v in zip(coef, ker):
    for m, cc in v.items(): E[m] = E.get(m, F(0)) + c_ * cc
E = {m: c for m, c in E.items() if c}; print("extra invariant: nonzero monomials", len(E), "; levels", sorted({excess(m) for m in E}), flush=True)
N = sp.Rational(-143, 4); t0 = time.time()
e = G.vev_polynomial(E, N); print(f"VEV computed in {time.time()-t0:.0f}s", flush=True)
P, Q = sp.symbols('P Q'); beta = 1 / sp.sqrt(q); Lam = lambda p: p**2 - Q**2 + sp.Rational((q - 1)**2, 2 * q)
pred = sp.expand(-sp.prod([Lam(P + sp.I * sp.sqrt(2) * m * beta) for m in range(-(q - 2) // 2, (q - 2) // 2 + 1)]))
print("eigenvalue (article momenta):", sp.factor(e))
print("prediction:", sp.factor(pred))
print("FUSION PREDICTION AT q = 8:", "CONFIRMED" if sp.expand(e - pred) == 0 else "REFUTED", "; difference:", sp.factor(sp.expand(e - pred)))
