"""Lead's blind comparison for EXC3: Fable's registered PREDICTION_EXC3.json (sha256 7991f7af..., snapshotted 03:51:27Z) vs the lead's own
Solution 2 level-one I_5 blocks (d1_data_exc3.py with5, run by the lead after registration).  Normalisation: matrix / [P^6] of the vacuum I_5
eigenvalue.  Usage: lead_compare_exc3.py [tamper]"""
import json, sys, hashlib, sympy as sp
P, Q = sp.symbols('P Q')
tamper = len(sys.argv) > 1
pred_path = 'PREDICTION_EXC3_snapshot.json'  # KIT PATCH: the lead's snapshot next to this script (was ../prereg_exc3/PREDICTION_EXC3.json, outside the archive); same sha256 7991f7af8dc89cc9...
assert hashlib.sha256(open(pred_path,'rb').read()).hexdigest().startswith('7991f7af8dc89cc9')
pred = json.load(open(pred_path))
loc = {'P': P, 'Q': Q}
bad = 0
for tf, tk in (('2', '2'), ('9_4', '9/4')):
    d = json.load(open(f'd1_sol2_t{tf}_with5.json'))
    b = d['blocks']['I5']
    vac = sp.expand(sp.sympify(b['vacuum'], locals=loc))
    c6 = sp.Poly(vac, P, Q).coeff_monomial(P**6)
    M = sp.Matrix([[sp.sympify(x, locals=loc) for x in row] for row in b['matrix']]) / c6
    if tamper: M[0, 0] += sp.Rational(1, 1000)
    tr, de = sp.expand(M.trace()), sp.expand(M.det())
    pf = pred['fibres'][tk if tk in pred['fibres'] else tf]
    ptr = sp.expand(sp.sympify(pf['I5_trace'], locals=loc)); pde = sp.expand(sp.sympify(pf['I5_det'], locals=loc))
    mt = 0 if sp.expand(tr - ptr) == 0 else len(sp.Poly(sp.expand(tr - ptr), P, Q).terms())
    md = 0 if sp.expand(de - pde) == 0 else len(sp.Poly(sp.expand(de - pde), P, Q).terms())
    nt, nd = len(sp.Poly(ptr, P, Q).terms()), len(sp.Poly(pde, P, Q).terms())
    print(f't = {tk}: trace mismatches {mt}/{nt}, det mismatches {md}/{nd}', flush=True)
    bad += mt + md
print('RESULT', 'MATCH' if bad == 0 else 'MISMATCH'); sys.exit(0 if bad == 0 else 2)
