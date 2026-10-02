"""SECOND2 main run v2: Suzuki class test for Sol 2 and Sol 3.
Spin 3 data NOT available for Sols 2/3 (only W4 for Sol 1). Fit on SPIN 5 instead.
H1 (Sol 2): free dict {sX, sY, cX, cY}, fit on spin 5 (4 monomials: X^2Y, X^2, XY, X).
H2a (Sol 3): same as H1.
H2b (Sol 3 rotated): A = sX*(X+pi)^2 + cX, 5 unknowns, fit on spin 5 (4 eqs) + spin 7 (1 eq).
Fibres: t = 2, 9/4, 4. M = (t+3)/(t-1) per fibre.
Report EVERY exact solution set per fibre BEFORE cross-fibre verdict.
Exit 0 always (report-only)."""
import sys, json; sys.dont_write_bytecode = True
from fractions import Fraction as F
import sympy as sp
from suzuki_wkb import riccati, charge

X, Y, t = sp.symbols('X Y t')
A, L = sp.symbols('A L')
sX, sY, cX, cY = sp.symbols('sX sY cX cY')
pi_rot = sp.Symbol('pi_rot')

def load_cert(sol, k, tv):
    """Certified top-layer for solution `sol`, weight W=2k, fibre tv."""
    TAB = {W: json.load(open(f'vev_sol{sol}_w{W}.json'))['coefficients'] for W in (6, 8, 10)}
    e = 0
    for kk, v in TAB[2*k].items():
        i, j = [int(z) for z in kk.replace('P^', '').replace('Q^', '').split()]
        e += sp.sympify(v, locals={'t': t}).subs(t, sp.Rational(tv.numerator, tv.denominator)) * X**(i//2) * Y**(j//2)
    e = sp.expand(e)
    lead = sp.Poly(e, X, Y).coeff_monomial(X**k)
    if lead == 0: return None
    return sp.expand(e / lead)

def wkb_poly(M_val, k):
    """WKB charge at R[2k] -> polynomial in A, L (unmapped)."""
    M = F(M_val)
    R = riccati(M, 10)
    cl = None
    for cls, (b0, f0, poly) in charge(R[2*k], M).items():
        if cls[0] != 'INTEGER-f':
            cl = (b0, f0, poly); break
    if cl is None:
        raise RuntimeError(f"No non-integer-f class at R[{2*k}], M={M}")
    b0, f0, poly = cl
    expr = sp.sympify(0)
    for (e, f), c in poly.items():
        e_int = int(e)
        assert e_int % 2 == 0, f"Odd e={e} at R[{2*k}], M={M}"
        expr += sp.Rational(c.numerator, c.denominator) * A**(e_int//2) * L**f
    return sp.expand(expr), (b0, f0)

def map_and_monic(expr_AL, k, rotated=False):
    """Map dictionary and monic-normalise. Returns mapped expression or None."""
    if rotated:
        mapped = sp.expand(expr_AL.subs({A: sX*(X + pi_rot)**2 + cX, L: sY*Y + cY}))
    else:
        mapped = sp.expand(expr_AL.subs({A: sX*X + cX, L: sY*Y + cY}))
    lead = sp.Poly(mapped, X, Y).coeff_monomial(X**k)
    if lead == 0: return None
    return sp.expand(mapped / lead)

def test_fibre(sol, tv, M_val, label, rotated=False):
    """Test one fibre. Returns list of (solution_dict, all_match, spin_details)."""
    print(f"\n--- {label}: Sol {sol}, t={tv}, M={M_val} ---", flush=True)

    # Fit on spin 5 (k=3): monomials X^2*Y, X^2, X*Y, X
    C5 = load_cert(sol, 3, tv)
    if C5 is None:
        print("  No certified spin-5 data; skipping.", flush=True)
        return []

    expr5, info5 = wkb_poly(M_val, 3)
    mapped5 = map_and_monic(expr5, 3, rotated)
    if mapped5 is None:
        print("  WKB spin-5 mapping failed.", flush=True)
        return []

    diff5 = sp.Poly(sp.expand(mapped5 - C5), X, Y)
    fit_monomials = (X**2 * Y, X**2, X * Y, X)
    fit_eqs = [sp.numer(sp.together(diff5.coeff_monomial(m))) for m in fit_monomials]

    if rotated:
        # 5 unknowns need 5 equations: add X^3*Y from spin 7
        C7 = load_cert(sol, 4, tv)
        if C7 is None:
            print("  No certified spin-7 data for 5th equation.", flush=True)
            return []
        expr7, _ = wkb_poly(M_val, 4)
        mapped7 = map_and_monic(expr7, 4, rotated)
        if mapped7 is None:
            print("  WKB spin-7 mapping failed.", flush=True)
            return []
        diff7 = sp.Poly(sp.expand(mapped7 - C7), X, Y)
        fit_eqs.append(sp.numer(sp.together(diff7.coeff_monomial(X**3 * Y))))
        unknowns = [sX, sY, cX, cY, pi_rot]
    else:
        unknowns = [sX, sY, cX, cY]

    print(f"  Fit on spin 5 (+spin 7 for H2b): {len(fit_eqs)} equations in {len(unknowns)} unknowns", flush=True)
    for i, eq in enumerate(fit_eqs):
        print(f"    eq[{i}]: {eq}", flush=True)

    try:
        sols = [so for so in sp.solve(fit_eqs, unknowns, dict=True) if so.get(sX, 1) != 0]
    except Exception as ex:
        print(f"  SOLVE FAILED: {type(ex).__name__}: {ex}", flush=True)
        return [("SOLVE_FAILED", None, {str(ex): True})]

    print(f"  Solutions found: {len(sols)}", flush=True)
    results = []
    for i, so in enumerate(sols):
        print(f"\n  Solution [{i}]:", flush=True)
        for u in unknowns:
            val = so[u]
            print(f"    {u.name} = {val}", flush=True)

        # Check ALL spins (5, 7, 9 — spin 3 not available)
        all_match = True
        spin_details = {}
        for k in (3, 4, 5):  # spins 5, 7, 9
            Ck = load_cert(sol, k, tv)
            if Ck is None:
                spin_details[2*k-1] = 'no data'; continue
            expr_k, _ = wkb_poly(M_val, k)
            mapped_k = map_and_monic(expr_k, k, rotated)
            if mapped_k is None:
                spin_details[2*k-1] = 'WKB fail'; all_match = False; continue
            mapped_k_sub = sp.expand(mapped_k.subs(so))
            lead_k = sp.Poly(mapped_k_sub, X, Y).coeff_monomial(X**k)
            if lead_k == 0:
                spin_details[2*k-1] = 'lead=0'; all_match = False; continue
            mapped_k_monic = sp.expand(mapped_k_sub / lead_k)
            diff_k = sp.Poly(sp.expand(mapped_k_monic - Ck), X, Y)
            nz = [(str(m), str(c)) for m, c in zip(diff_k.monoms(), diff_k.coeffs()) if c != 0]
            spin_details[2*k-1] = 'MATCH' if not nz else f'{len(nz)} MISMATCH'
            if nz:
                all_match = False
                print(f"    spin {2*k-1}: {nz[0][0]} -> {nz[0][1][:80]}...", flush=True)
            else:
                print(f"    spin {2*k-1}: MATCH", flush=True)

        results.append((so, all_match, spin_details))
        status = 'ALL MATCH' if all_match else 'FAILS'
        print(f"  => Solution [{i}]: {status}", flush=True)

    return results

def main():
    print("=" * 70)
    print("SECOND2: Suzuki class test for Sol 2 and Sol 3")
    print("Fit on spin 5 (spin 3 data unavailable for Sols 2/3)")
    print("=" * 70, flush=True)

    fibres = [F(2), F(9, 4), F(4)]
    M_vals = [(tv + 3) / (tv - 1) for tv in fibres]
    for tv, Mv in zip(fibres, M_vals):
        print(f"Fibre t={tv}: M={Mv}")

    sol2_results = {}
    sol3a_results = {}
    sol3b_results = {}

    # === H1: SOL 2 ===
    print("\n" + "=" * 70)
    print("H1: SOL 2 (standard dictionary A = sX*X + cX)")
    print("=" * 70, flush=True)
    for tv, Mv in zip(fibres, M_vals):
        try:
            res = test_fibre(2, tv, Mv, "H1 Sol2", rotated=False)
            sol2_results[str(tv)] = res
        except ZeroDivisionError as ex:
            print(f"  DEGENERATE at t={tv}: {ex}", flush=True)
            sol2_results[str(tv)] = [("DEGENERATE", None, str(ex))]
        except Exception as ex:
            import traceback; traceback.print_exc()
            sol2_results[str(tv)] = [("ERROR", None, f"{type(ex).__name__}: {ex}")]

    # === H2a: SOL 3 standard ===
    print("\n" + "=" * 70)
    print("H2a: SOL 3 (standard dictionary A = sX*X + cX)")
    print("=" * 70, flush=True)
    for tv, Mv in zip(fibres, M_vals):
        try:
            res = test_fibre(3, tv, Mv, "H2a Sol3", rotated=False)
            sol3a_results[str(tv)] = res
        except ZeroDivisionError as ex:
            print(f"  DEGENERATE at t={tv}: {ex}", flush=True)
            sol3a_results[str(tv)] = [("DEGENERATE", None, str(ex))]
        except Exception as ex:
            import traceback; traceback.print_exc()
            sol3a_results[str(tv)] = [("ERROR", None, f"{type(ex).__name__}: {ex}")]

    # === H2b: SOL 3 rotated ===
    print("\n" + "=" * 70)
    print("H2b: SOL 3 (ROTATED dictionary A = sX*(X+pi)^2 + cX)")
    print("=" * 70, flush=True)
    for tv, Mv in zip(fibres, M_vals):
        try:
            res = test_fibre(3, tv, Mv, "H2b Sol3 rot", rotated=True)
            sol3b_results[str(tv)] = res
        except ZeroDivisionError as ex:
            print(f"  DEGENERATE at t={tv}: {ex}", flush=True)
            sol3b_results[str(tv)] = [("DEGENERATE", None, str(ex))]
        except Exception as ex:
            import traceback; traceback.print_exc()
            sol3b_results[str(tv)] = [("ERROR", None, f"{type(ex).__name__}: {ex}")]

    # === CROSS-FIBRE CONSISTENCY VERDICT ===
    print("\n" + "=" * 70)
    print("CROSS-FIBRE CONSISTENCY VERDICT")
    print("=" * 70, flush=True)

    for name, results in [("SOL 2 (H1)", sol2_results),
                           ("SOL 3 standard (H2a)", sol3a_results),
                           ("SOL 3 rotated (H2b)", sol3b_results)]:
        matching_fibres = []
        all_solutions_reported = []
        for tv_str, res in sorted(results.items()):
            n_sols = len(res)
            matches = [r for r in res if isinstance(r, tuple) and len(r) >= 2 and r[1] is True]
            if matches:
                matching_fibres.append(tv_str)
            # Report all solutions
            for i, item in enumerate(res):
                if isinstance(item, tuple) and len(item) >= 2 and isinstance(item[0], dict):
                    so, all_m, details = item
                    sol_desc = {u.name: str(so.get(u, '?')) for u in ([sX,sY,cX,cY,pi_rot] if 'H2b' in name else [sX,sY,cX,cY]) if u in so}
                    all_solutions_reported.append(f"    t={tv_str} sol[{i}]: {sol_desc} -> {'MATCH' if all_m else 'FAIL'}")
        print(f"\n{name}:")
        for line in all_solutions_reported:
            print(line, flush=True)
        if len(matching_fibres) == len(fibres):
            verdict = "PASS: full match at ALL fibres"
        elif matching_fibres:
            verdict = f"PARTIAL: matches only at {matching_fibres}"
        else:
            verdict = "FAIL: no fibre has a full match"
        print(f"  VERDICT: {verdict}", flush=True)

    print("\ndone", flush=True)

if __name__ == '__main__':
    main()
