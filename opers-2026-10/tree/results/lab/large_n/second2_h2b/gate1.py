"""SECOND2 GATE 1 v4: Engine validation.
Convention: WKB poly in A,L -> map dict A=sX*X+cX, L=sY*Y+cY -> monic in X,Y (closed_dict.py convention).
M fixed to known value (5 for plant; 51/10 for tamper). Fit on spin 3, check 5,7,9.
PASS criterion: AT LEAST ONE solution from the fit matches ALL higher-spin checks.
Exit 0 iff plant PASSES (some solution matches all) AND tamper FAILS (no solution matches all)."""
import sys, json; sys.dont_write_bytecode = True
from fractions import Fraction as F
import sympy as sp
from suzuki_wkb import riccati, charge

X, Y, t = sp.symbols('X Y t')
A, L = sp.symbols('A L')
sX, sY, cX, cY = sp.symbols('sX sY cX cY')

W4 = json.load(open('vev_w4_points.json'))['1']
TAB = {W: json.load(open(f'vev_sol1_w{W}.json'))['coefficients'] for W in (6, 8, 10)}

def cert(k, tv):
    if k == 2:
        d = W4.get(str(F(tv.numerator, tv.denominator)))
        if d is None: return None
        e = sum(sp.Rational(v) * X**int(kk.split(',')[0]) * Y**int(kk.split(',')[1]) for kk, v in d.items())
    else:
        e = 0
        for kk, v in TAB[2*k].items():
            i, j = [int(z) for z in kk.replace('P^', '').replace('Q^', '').split()]
            e += sp.sympify(v, locals={'t': t}).subs(t, sp.Rational(tv.numerator, tv.denominator)) * X**(i//2) * Y**(j//2)
    e = sp.expand(e)
    lead = sp.Poly(e, X, Y).coeff_monomial(X**k)
    return sp.expand(e / lead)

def wkb_mapped(M_val, k):
    M = F(M_val)
    R = riccati(M, 10)
    cl = None
    for cls, (b0, f0, poly) in charge(R[2*k], M).items():
        if cls[0] != 'INTEGER-f':
            cl = (b0, f0, poly); break
    if cl is None:
        raise RuntimeError(f"No non-integer-f class at R[{2*k}], M={M}")
    b0, f0, poly = cl
    expr_AL = sp.sympify(0)
    for (e, f), c in poly.items():
        e_int = int(e)
        assert e_int % 2 == 0, f"Odd e={e} at R[{2*k}], M={M}"
        expr_AL += sp.Rational(c.numerator, c.denominator) * A**(e_int//2) * L**f
    expr_AL = sp.expand(expr_AL)
    mapped = sp.expand(expr_AL.subs({A: sX*X + cX, L: sY*Y + cY}))
    lead = sp.Poly(mapped, X, Y).coeff_monomial(X**k)
    return sp.expand(mapped / lead)

def check_solution(so, M_val, tv):
    """Check if a dictionary solution matches all spins. Returns (all_match, details)."""
    all_match = True
    details = []
    for k in (3, 4, 5):
        Ck = cert(k, tv)
        if Ck is None: continue
        mapped_k = wkb_mapped(M_val, k).subs(so)
        lead_k = sp.Poly(sp.expand(mapped_k), X, Y).coeff_monomial(X**k)
        mapped_k = sp.expand(mapped_k / lead_k)
        diff_k = sp.Poly(sp.expand(mapped_k - Ck), X, Y)
        nz = [(str(m), str(c)) for m, c in zip(diff_k.monoms(), diff_k.coeffs()) if c != 0]
        status = 'MATCH' if not nz else f'{len(nz)} MISMATCH'
        if nz: all_match = False
        details.append(f"spin {2*k-1}: {status}")
    return all_match, details

def run_case(M_val, label):
    tv = F(2)
    print(f"\n{'='*60}")
    print(f"GATE 1 {label}: M = {M_val}, t = 2")
    print(f"{'='*60}", flush=True)

    C3 = cert(2, tv)
    mapped3 = wkb_mapped(M_val, 2)
    diff3 = sp.Poly(sp.expand(mapped3 - C3), X, Y)
    fit_eqs = [sp.numer(sp.together(diff3.coeff_monomial(m))) for m in (X*Y, X, Y, 1)]
    print(f"Spin 3 fit equations:", flush=True)
    for eq in fit_eqs: print(f"  {eq}", flush=True)

    try:
        sols = [so for so in sp.solve(fit_eqs, [sX, sY, cX, cY], dict=True) if so.get(sX, 1) != 0]
    except Exception as ex:
        print(f"  SOLVE FAILED: {ex}", flush=True)
        return False
    print(f"  Solutions: {len(sols)}", flush=True)
    any_pass = False
    for i, so in enumerate(sols):
        print(f"\n  Solution [{i}]: {so}", flush=True)
        all_match, details = check_solution(so, M_val, tv)
        for d in details: print(f"    {d}", flush=True)
        print(f"    => {'MATCH ALL' if all_match else 'FAILS'}", flush=True)
        if all_match:
            any_pass = True
            # Report the dictionary in closed form for comparison with SUPER3b
            print(f"    Dictionary: sX={so[sX]}, sY={so[sY]}, cX={so[cX]}, cY={so[cY]}")
            print(f"    (SUPER3b expects: sX=-4M(M+1)={-4*M_val*(M_val+1)}, sY=-(M+1)={-(M_val+1)}, cX=0, cY=(M-1)^2/4={(M_val-1)**2/4})")

    result = 'PASS' if any_pass else 'FAIL'
    print(f"\n  RESULT: {result}", flush=True)
    return any_pass

plant_ok = run_case(F(5), "PLANT")
tamper_ok = run_case(F(51, 10), "TAMPER")

print(f"\n{'='*60}")
if plant_ok and not tamper_ok:
    print("GATE 1: PASS (plant recovered, tamper fails as expected)")
    sys.exit(0)
else:
    print(f"GATE 1: FAIL (plant={'ok' if plant_ok else 'fail'}, tamper={'ok(unexpected)' if tamper_ok else 'fail(correct)'})")
    sys.exit(1)
