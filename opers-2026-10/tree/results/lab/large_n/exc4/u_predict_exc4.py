"""EXC4 B3: assemble PREDICTION_EXC4.json from the four-point runs at the hold-out points (t = 9/4 and further points): for each point, the monic
I_3 eliminant over the sigma-invariant four-point solutions (predicted = the quintic factor of the level-two singlet I_3 block's characteristic
polynomial, the sixth factor being the level-shift state e = vacuum at Q^2 -> Q^2 + 4) and the I_1 eliminant.  Opens NO level-two data at
those points.  Usage: u_predict_exc4.py <t> <P> <Q> [<t> <P> <Q> ...]"""
import json, re, sys
import sympy as sp
sys.dont_write_bytecode = True
x = sp.Symbol('x')
out = {'stage': 'EXC4 part B', 'seat': 'Fable (secondary)', 'object': 'Solution 3, level-two singlet block (6 states) of I_3: the degree-5 factor of the characteristic polynomial of the block (eigenvalues divided by the P^4 coefficient of the vacuum I_3), the sixth root being the level-shift state (vacuum at Q^2 -> Q^2 + 4)',
       'variables': 'bare P = P_X, Q^2 = pi^2 - rho^2 as in lab/excited_states.py; points given as rational (P, Q)', 'source': 'four apparent singularities +-z1, +-z3 of type {-1,1,2,4}, sigma-invariant pairs; I_3 from the effective quartic', 'points': []}
args = sys.argv[1:]
for i in range(0, len(args), 3):
    tstr, Ps, Qs = args[i:i + 3]
    TAG = 't%s_P%s_Q%s' % (tstr.replace('/', '_'), Ps.replace('/', '-'), Qs.replace('/', '-'))
    res = open('h_four_%s_four.out' % TAG).read()
    vd = int(re.search(r'vdim:\n(-?\d+)', res).group(1))
    def mon(key, nm):
        p_ = sp.Poly(sp.sympify(re.search(r'%s\[1\]=(.*)' % key, res).group(1).replace('^', '**').replace(nm, 'x'), locals={'x': x}), x)
        return sp.Poly(p_.as_expr() / p_.LC(), x)
    pe_, pi_ = mon('Je', 'ee'), mon('Ji', 'ii')
    out['points'].append({'t': tstr, 'P': Ps, 'Q': Qs, 'vdim_ordered': vd, 'I3_quintic_monic': str(pe_.as_expr()).replace('x', 'e'), 'I1_eliminant': str(pi_.as_expr()).replace('x', 'e')})
    print(TAG, 'vdim', vd, 'I_3 degree', pe_.degree())
json.dump(out, open('PREDICTION_EXC4.json', 'w'), indent=1)
print('written PREDICTION_EXC4.json with %d points' % len(out['points']))
