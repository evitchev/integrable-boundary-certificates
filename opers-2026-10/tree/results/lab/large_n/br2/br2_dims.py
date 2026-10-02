"""BR2 (4): boundary dimensions Delta = h(v) = v.v/2 - v.Q (Q = (0, rho), rho = p/2 - 1/p; doubled chiral weight = boundary dimension, BLZ (5.1) fn. 8 convention),
for B_+- of each solution (and Sol 3's dual pair) and B_vir; relevance ranges 0 < Delta < 1 in t (p^2 = 2(t-1)/(t+1)).  Also the LVZ hairpin parameter n of each pair."""
import sys; sys.dont_write_bytecode = True
import sympy as sp
t = sp.Symbol('t', real=True); P2 = sp.Symbol('P2')     # P2 = p^2
p = sp.sqrt(P2); rho = p / 2 - 1 / p
data = {'Sol 1': (1 - P2 / 4, -p / 2), 'Sol 2': (1 - P2, -p), 'Sol 3': (1 / P2, 1 / p), "Sol 3 dual": (P2, -p), 'Virasoro e^(p phi)': (0, p)}
p2t = 2 * (t - 1) / (t + 1)
for lab, (a2, beta) in data.items():
    Delta = sp.simplify(a2 / 2 + beta**2 / 2 - beta * rho)
    Dt = sp.factor(sp.simplify(Delta.subs(P2, p2t)))
    rel = sp.reduce_inequalities([Dt > 0, Dt < 1], t) if Dt.free_symbols else ('marginal' if Dt == 1 else Dt)
    norm = sp.simplify(a2 + beta**2); cross = sp.simplify(-a2 + beta**2)
    hair = sp.simplify(-(cross + 1)) if sp.simplify(norm - 1) == 0 else 'not norm 1 (no LVZ hairpin match)'
    print(f'{lab}: Delta = {sp.factor(Delta)} = {Dt};  relevant (0 < Delta < 1) for: {rel};  pair norm {norm}, LVZ hairpin n = -(v+.v- + 1) = {hair}')
