"""PILLOW1 route (3): Newton polygon at p = infinity of L4 = ((p^2 + E) chi)'''' - lam chi'' + alpha chi' - chi (Laplace of Suzuki(M = -2) x^4).
Route (1d): exponential gauges exp(beta s^h), h = 1, 2, 3, on the pillow's 3-term and 5-term forms: can any beta reduce or change the shift pattern to Sol 1's / Sol 2's?"""
import sys; sys.dont_write_bytecode = True
import sympy as sp
p, E, lam, al, q = sp.symbols('p E lam alpha q')
chi = sp.Function('chi')
# Newton polygon: write each term as (derivative order j, power of p after expressing d/dp = p^-1 theta): p^m d^j -> weight m - j
terms = [(4, 2), (4, 0), (3, 1), (2, 0), (2, 0), (1, 0), (0, 0)]   # from expanding ((p^2+E) chi)'''' = (p^2+E) chi'''' + 8 p chi''' + 12 chi'' ; - lam chi''; alpha chi'; -chi
pts = sorted({(j, m - j) for j, m in terms})
print('Newton points (order j, p-weight m - j):', pts)
# at p -> infinity the polygon is the UPPER hull of the weights; slope-0 (regular) solutions need the top weight attained at more than one order j
top = max(w for _, w in pts); at_top = [j for j, w in pts if w == top]
print('maximal p-weight', top, 'attained at orders', at_top, '-> regular-type (slope 0) formal solutions at infinity exist only if the maximum is attained at order >= 1: ',
      max(at_top) >= 1)
print('edge from (0, 0) to (4, -2) has slope 1/2 and contains (1, -1) and (3, -2)? -> all four formal solutions are exp(c p^(1/2)) type (rank 1/2): NO slope-0 part')
# Route (1d): exponential gauges on the pillow 3-term form
s, X, Y, kap, beta, T = sp.symbols('s X Y kappa beta T')
sq = sp.sqrt(2) * sp.sqrt(-X)
p0 = lambda TT: (2 * TT - 1) * (2 * TT + 1) / 4
p2 = lambda TT: -((2 * TT + 2 * sq + 2)**2 - (4 * Y + 9)) / 4
def op_terms(h, form):
    # operator as dict shift -> polynomial in T (theta to the right); exponential gauge: theta -> theta + beta h s^h, applied as conjugation
    y = sp.Function('y')
    if form == 3:
        L = lambda f: p0(0) * f + 0
    G = sp.exp(beta * s**h)
    def theta(f): return s * sp.diff(f, s)
    def Lop(f):
        if form == 3:
            c = {0: p0, 2: p2}
            out = sum(s**m * (sp.Poly(c[m](T), T).coeff_monomial(T**2) * theta(theta(f)) + sp.Poly(c[m](T), T).coeff_monomial(T) * theta(f) + sp.Poly(c[m](T), T).coeff_monomial(1) * f) for m in c)
            return out - 4 * kap**2 * s**3 * f
    expr = sp.expand(sp.simplify(Lop(G * y(s)) / G))
    d2 = sp.Derivative(y(s), (s, 2)); d1 = sp.Derivative(y(s), s)
    a2 = expr.coeff(d2); a1 = sp.expand(expr - a2 * d2).coeff(d1); a0 = sp.expand(expr - a2 * d2 - a1 * d1).coeff(y(s))
    mons = set()
    for cc in (sp.expand(a2 * s**2), sp.expand(a1 * s), sp.expand(a0)):
        mons |= {m[0] for m in sp.Poly(cc, s).monoms()} if cc != 0 else set()
    return sorted(mons)
for h in (1, 2, 3):
    print(f'exp(beta s^{h}) gauge on the pillow 3-term form: shifts present (generic beta): {op_terms(h, 3)}')
