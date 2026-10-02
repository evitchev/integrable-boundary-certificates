"""EXC1 O3 (hold-out): prediction of the characteristic polynomial of I_5 on the level-1 singlet block of Sol 1 from the oper
with ONE apparent singularity (trivial-monodromy quadratic), at t = 2 and t = 9/4.  Opens no excited-state data.
Normalisation: eigenvalue divided by the PX^6 coefficient of the vacuum eigenvalue; variables PX, pi.
-> PREDICTION_EXC1.json"""
import hashlib, json, sys
import sympy as sp
sys.dont_write_bytecode = True
assert hashlib.sha256(open('SEAL_EXC1.md', 'rb').read()).hexdigest() == open('SEAL_EXC1.sha256').read().split()[0]
l0, l1, z1, e, PX, PI = sp.symbols('l0 l1 z1 e PX PI')
out = {'seal': open('SEAL_EXC1.sha256').read().split()[0], 'object': 'monic characteristic polynomial in e of I_5/(PX^6 vacuum coefficient) on the level-1 singlet block, Sol 1', 'cells': {}}
for tstr in ('2', '9/4'):
    t = sp.Rational(tstr); M = (t + 3) / (t - 1)
    od = json.load(open('o1_t%s_n1.json' % tstr.replace('/', '_')))
    R6 = sp.sympify(od['R']['6'], locals={'Symbol': sp.Symbol})
    c6 = sp.Poly(R6, l0, l1, z1).coeff_monomial(l1**6)
    vac_lead = sp.Poly(sp.expand((R6 / c6).subs(z1, 0)), l0, l1).coeff_monomial(l0**6) * (4 * M * (M + 1))**3
    L5 = sp.expand(R6 / c6 / vac_lead)
    cond = 2 * M * z1**2 + l0 * (M - 1) * z1 - 2 * l1**2 + (M - 1)**2 / 2
    res = sp.expand(sp.resultant(sp.expand(cond), sp.expand(e - L5), z1))
    resP = sp.Poly(res, l0)
    assert all(k[0] % 2 == 0 for k in resP.as_dict())
    res_m = sum(c * (4 * M * (M + 1) * PX**2) ** (k[0] // 2) for k, c in resP.as_dict().items())
    res_m = sp.Poly(sp.expand(res_m.subs(l1**2, (M + 1) * PI**2)), e)
    res_m = sp.Poly(sp.expand(res_m.as_expr() / res_m.LC()), e)
    out['cells']['t=' + tstr] = {'M': str(M), 'charpoly': str(res_m.as_expr()), 'coeff_e1': str(res_m.coeff_monomial(e)), 'coeff_e0': str(res_m.coeff_monomial(1))}
    print('t = %s: predicted I_5 level-1 characteristic polynomial: degree 2 in e; trace = %s' % (tstr, sp.factor(-res_m.coeff_monomial(e))))
json.dump(out, open('PREDICTION_EXC1.json', 'w'), indent=1, sort_keys=True)
print(hashlib.sha256(open('PREDICTION_EXC1.json', 'rb').read()).hexdigest(), ' PREDICTION_EXC1.json')
