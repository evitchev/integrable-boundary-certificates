"""WUW1 post-hoc (labelled): the discriminating version of P2.  Weight-3 fields in one boson: basis f1 = (dU)^3, f2 = dU d^2U, f3 = d^3U (same for W).
Find all 3x3 M with sum M_ij f_i(U) g_j(W) in K_6 + d (linear), and the rank-1 members (products).  Is W_U W_W the unique product (up to scale)?
Also: tampered product S_3(21/10)_U S_3(21/10)_W in K_6?"""
import os, sys
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'wuw1.py')).read()
pre = src[:src.index('# ================= (1) the weight-3 objects')].replace("open(os.path.join(HERE, 'run_wuw1.log'), 'w')", "open(os.path.join(HERE, 'run_posthoc_products.log'), 'w')")
exec(pre)
SCR = [{U: F(1, 2), Wc: F(1, 2)}, {U: F(2), Wc: F(0)}, {U: F(0), Wc: F(2)}]
I3 = reduce_mod_d(commutant(4, SCR), 4)[0]; K6 = ad_kernel(I3, 6)
fU = [{((U, 0),) * 3: F(1)}, {((U, 0), (U, 1)): F(1)}, {((U, 2),): F(1)}]
fW = [{((Wc, 0),) * 3: F(1)}, {((Wc, 0), (Wc, 1)): F(1)}, {((Wc, 2),): F(1)}]
prods = [mul(fU[i], fW[j]) for i in range(3) for j in range(3)]
dims = [dpoly({m: F(1)}) for m in monomials(5)]
keys = sorted({k for f in prods + K6 + dims for k in f}); idx = {k: n for n, k in enumerate(keys)}
cols = []
for f in prods + K6 + dims: cols.append({idx[k]: v for k, v in f.items()})
ns = nullspace(cols, len(keys))
Ms = [x[:9] for x in ns if any(v != 0 for v in x[:9])]
Mb, _ = rref([list(m) for m in Ms], 9); say(f'linear space of M (3x3, f_i(U) g_j(W)) with product in K_6 + d: dim {len(Mb)}')
for r in Mb: say(f'   basis row (M_11..M_33): {[str(v) for v in r]}')
xs = sp.symbols('x0:%d' % len(Mb)); M = sp.Matrix(3, 3, lambda i, j: sum(xs[k] * sp.Rational(Mb[k][3 * i + j].numerator, Mb[k][3 * i + j].denominator) for k in range(len(Mb))))
minors = [sp.expand(M.extract([a, b], [c, d]).det()) for a in range(3) for b in range(a + 1, 3) for c in range(3) for d in range(c + 1, 3)]
sols = sp.solve([m for m in minors if m != 0], xs, dict=True) if len(Mb) > 1 else [{}]
say(f'rank-1 (product) solutions: {sols}')
WUv = [F(4, 3), F(2), F(1, 3)]   # W_U coefficients on (f1, f2, f3)
target = [WUv[i] * WUv[j] for i in range(3) for j in range(3)]
say(f'W_U W_W as M: {[str(v) for v in target]}')
for so in sols:
    Mv = M.subs(so); say(f'   solution M = {Mv.tolist()}; rank {Mv.rank()}')
tg = mul(schur(F(21, 10), U, 3), schur(F(21, 10), Wc, 3)); say(f'TAMPER S_3(21/10)_U S_3(21/10)_W in K_6 + d: {in_span_mod_d(tg, K6, 6)} (must be False)')
say(f'(dU)^3 (dW)^3 in K_6 + d: {in_span_mod_d(prods[0], K6, 6)}')
