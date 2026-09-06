"""First-order certificate (2026-09-04): around the N = 1 anchor the cylindrical
family's deformation is NOT produced by any scalar Schroedinger-type
deformation of the paperclip operator in its own variable, while the
pipeline that says so reproduces the paperclip's own deformation exactly.

Setting.  At t -> infinity (N = 1) the cylindrical VEVs of Solutions 1 and 3
equal the paperclip's (hep-th/0404195 sec. 4) at n = -1 up to
normalization, and the paperclip ODE (62) at n = -1,
    Psi'' = [kappa^2 (1+u)^m + V0] Psi,  V0 = -X/2 u/(1+u) + (Y/2 + 1/4) u/(1+u)^2,
X = P^2, Y = Q^2 (article momenta), m = -1, reproduces them through its
WKB coefficients S_s (sec. 5 prescription, code/wkb.py).  Expanding in
e = 1/t, the first-order change of S_s under a deformation
    delta m,  delta V = sum_k c_k(X, Y) f_k(u)  (c_k affine),
plus the free normalization changes delta alpha_s I_s + delta beta_s, is
LINEAR in the unknowns; the response column of each term f_k is
d S_s/d c_k at the base.  The columns (spins 1, 3, 5, 7, 9) are PINNED
artifacts (results/lab/first_order/, rational polynomials in X, Y; the
log-2 constant of the origin-log class kept as the symbol TT), produced in
the lab with an independent Wolfram port of the WKB integrator; the
first-order data are the 1/t-expansion of the certified VEV tables
(weights 6, 8, 10: cyl_vev_tables) and of the weight-4 table pinned here
(vev_w4_points.json, the same recipe at 30 fibres, fitted exactly).

Checks:
  P1  independent recomputation: the pinned columns of delta m and of
      u/(1+u), u/(1+u)^2, u/(1+u)^3 at spins 1-5 equal those computed here
      with code/wkb.py (m symbolic, derivative at -1) -- exactly;
  P2  base at spins 3, 5, 7: the pinned base S_s^0 equals alpha_s I_s^pc + beta_s
      with the paper's alpha_s (validated in cyl_wkb_validate) and the
      cylindrical tables at t = infinity; at spin 9 the base is
      proportional to the table (alpha_9 = -1/1440): the paperclip
      correspondence at spin 9 (new);
  P3  POSITIVE CONTROL: with the paperclip's own n-deformation as data
      the system is consistent with exactly delta m = 1, delta A = X/2,
      delta B = Y/2 and no normalization change;
  P4  for Solutions 3 and 1, with spins 1-9 (55 rational data coordinates)
      and the spin-1 equation (delta I_1 = 0, the central charge moving
      only at second order), every listed class is INCONSISTENT:
      the two LVZ terms; the (1+u)-library u/(1+u)^b (b <= 6) with growth
      terms u, u^2; the library with the moving-exponent log columns at
      u = -1 and at the origin; and (spins 1-7, pinned columns) a fixed
      fourth point at u = -v0, v0 = 3, 9, 1/2, 2, with its transcendental
      as an indeterminate;
  P5  negative controls: the positive-control expectation with delta m
      replaced by 2 fails, and a tampered data coordinate makes the
      positive control inconsistent (the guards are live).
Fail-closed: any failure exits 1 without the CERTIFIED marker.
Runtime ~3-4 min (sympy rank computations).
"""
import sys, json, os, time
import sympy as sp
from wkb import wronskian_S, u, m as msym

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FO = os.path.join(ROOT, "results", "lab", "first_order")
fails = []
def require(cond, msg):
    print(("PASS  " if cond else "FAIL  ") + msg, flush=True)
    if not cond:
        fails.append(msg)

n, P, Q, t, e, X, Y, TT = sp.symbols('n P Q t e X Y TT')
# the paper's VEVs (sec. 4) and normalizations, as validated in cyl_wkb_validate.py
IPC = {3: (4+3*n)/(6*(2+n))*P**4 + (2+3*n)/(6*n)*Q**4 + P**2*Q**2 + (3+2*n)/(6*(2+n))*P**2 + (1+2*n)/(6*n)*Q**2 + (11+36*n+18*n**2)/(360*n*(2+n)),
       5: (-(6+5*n)*(8+5*n)/(120*(2+n)**2)*P**6 - (2+5*n)*(4+5*n)/(120*n**2)*Q**6 - (8+5*n)/(8*(2+n))*P**4*Q**2 - (2+5*n)/(8*n)*P**2*Q**4
           - (4+3*n)*(8+5*n)/(48*(2+n)**2)*P**4 - (2+3*n)*(2+5*n)/(48*n**2)*Q**4 - (2+10*n+5*n**2)/(8*n*(2+n))*P**2*Q**2
           - (30+225*n+250*n**2+76*n**3)/(480*n*(2+n)**2)*P**2 - (28+137*n+206*n**2+76*n**3)/(480*n**2*(2+n))*Q**2
           - (564+3410*n+7385*n**2+5680*n**3+1420*n**4)/(60480*n**2*(2+n)**2)),
       7: (n*(8+7*n)*(10+7*n)/(42*(2+n)**2*(2+7*n))*P**8 + (2+n)*(4+7*n)*(6+7*n)/(42*n**2*(12+7*n))*Q**8
           + 2*n*(10+7*n)/(3*(2+n)*(2+7*n))*P**6*Q**2 + 2*(2+n)*(4+7*n)/(3*n*(12+7*n))*P**2*Q**6 + P**4*Q**4
           + n*(5+4*n)*(10+7*n)/(9*(2+n)**2*(2+7*n))*P**6 + (2+n)*(3+4*n)*(4+7*n)/(9*n**2*(12+7*n))*Q**6
           + (6+49*n+28*n**2)/(3*(2+n)*(2+7*n))*P**4*Q**2 + (20+63*n+28*n**2)/(3*n*(12+7*n))*P**2*Q**4
           + (122+1407*n+1792*n**2+612*n**3)/(180*(2+n)**2*(2+7*n))*P**4 + (420+1583*n+1880*n**2+612*n**3)/(180*n**2*(12+7*n))*Q**4
           + 7*(40+338*n+985*n**2+816*n**3+204*n**4)/(30*n*(2+n)*(2+7*n)*(12+7*n))*P**2*Q**2
           + (2520+26166*n+102459*n**2+136612*n**3+74676*n**4+14552*n**5)/(1260*n*(2+n)**2*(2+7*n)*(12+7*n))*P**2
           + (3720+30202*n+89149*n**2+121284*n**3+70844*n**4+14552*n**5)/(1260*n**2*(2+n)*(2+7*n)*(12+7*n))*Q**2
           + (68760+632142*n+2264647*n**2+4095840*n**3+3708040*n**4+1610448*n**5+268408*n**6)/(151200*n**2*(2+n)**2*(2+7*n)*(12+7*n)))}
IPC = {k: sp.expand(v.subs({P: sp.sqrt(X), Q: sp.sqrt(Y)})) for k, v in IPC.items()}
ALPHA = {3: -n*(n+2)/(4*(3*n+2)*(3*n+4)), 5: -3*n**2*(n+2)**2/((5*n+2)*(5*n+4)*(5*n+6)*(5*n+8)), 7: -45*n**2*(n+2)**2/(32*(7*n+4)*(7*n+6)*(7*n+8)*(7*n+10))}
BETA = {3: sp.Rational(-1, 2880), 5: sp.Rational(1, 40320), 7: sp.Rational(-1, 215040)}
N0 = -1
LAM = {3: 6, 5: 40, 7: 70}      # the tables' normalization relative to the paper's I_s at n = -1 (weights 4, 6, 8)

def load_col(name):
    d = {}
    with open(os.path.join(FO, f"{name}.txt")) as fh:
        for line in fh:
            k, expr = line.split(" = ", 1)
            d[int(k)] = sp.sympify(eval(expr, sp.__dict__))
    return d

def fit_rational(samples, deg):
    ts = sorted(samples); rows = [[tv**k for k in range(deg + 1)] + [-samples[tv] * tv**k for k in range(deg + 1)] for tv in ts]
    ker = sp.Matrix(rows).nullspace()
    if not ker:
        return None
    v = ker[0]; f = sp.cancel(sum(v[k] * t**k for k in range(deg + 1)) / sum(v[deg + 1 + k] * t**k for k in range(deg + 1)))
    return f if all(f.subs(t, tv) == samples[tv] for tv in ts) else None

def cyl(sol, W):
    if W == 4:
        d = json.load(open(os.path.join(FO, "vev_w4_points.json")))[str(sol)]; expr = 0
        for key in sorted({k for v in d.values() if isinstance(v, dict) for k in v}):
            samples = {sp.Rational(tv): sp.Rational(v[key]) for tv, v in d.items() if isinstance(v, dict)}
            f = next(f_ for f_ in (fit_rational(samples, dg) for dg in range(0, 8)) if f_ is not None)
            a, b = [int(z) for z in key.split(',')]; expr += f * X**a * Y**b
        return expr
    d = json.load(open(os.path.join(ROOT, "results", "lab", "vev", f"vev_sol{sol}_w{W}.json"))); expr = 0
    for k, v in d['coefficients'].items():
        i, j = [int(z) for z in k.replace('P^', '').replace('Q^', '').split()]; expr += sp.sympify(v) * X**(i // 2) * Y**(j // 2)
    return expr

def first_order(sol, spins):
    """data[s] = (alpha_s^0, I_s^cyl(inf), dI_s^cyl/d(1/t))"""
    out = {}
    for s_ in spins:
        Ic = cyl(sol, s_ + 1); ser = sp.expand(sp.series(Ic.subs(t, 1 / e), e, 0, 2).removeO())
        I0 = sp.expand(ser.subs(e, 0)); I1 = sp.expand(sp.diff(ser, e).subs(e, 0))
        if s_ == 9:
            al = sp.Rational(json.load(open(os.path.join(FO, f"alpha9_sol{sol}.json")))['alpha9_0'])
        else:
            pc0 = sp.expand(IPC[s_].subs(n, N0)); lam0 = sp.simplify(sp.Poly(I0, X, Y).LC() / sp.Poly(pc0, X, Y).LC())
            if sp.expand(I0 - lam0 * pc0) != 0:
                raise RuntimeError(f"cylindrical table at t = infinity is not a multiple of the paperclip at spin {s_}")
            al = ALPHA[s_].subs(n, N0) / lam0
        out[s_] = (al, I0, I1)
    return out

def solve(cols, data, terms, spins, spin1=True):
    dm = sp.Symbol('dm'); unk = [dm]; expr = {s_: cols['m'][s_] * dm for s_ in (1,) + tuple(spins)}
    hasT = any(c_.has(TT) for nm in terms for c_ in cols[nm].values())
    for nm in terms:
        ks = sp.symbols(f'c0_{nm} c1_{nm} c2_{nm}'); unk += list(ks)
        for s_ in (1,) + tuple(spins): expr[s_] += cols[nm][s_] * (ks[0] + ks[1] * X + ks[2] * Y)
    eqs = []
    def add(dd):
        eqs.extend(sp.Poly(dd, X, Y, TT).coeffs() if hasT else sp.Poly(dd, X, Y).coeffs())
    if spin1:
        dal, dbe, dalT, dbeT = sp.symbols('dal1 dbe1 dal1T dbe1T'); unk += [dal, dbe] + ([dalT, dbeT] if hasT else [])
        add(sp.expand(expr[1] - ((dal + (TT * dalT if hasT else 0)) * (X + Y) + dbe + (TT * dbeT if hasT else 0))))
    for s_ in spins:
        al0, I0, I1 = data[s_]; dal, dbe, dalT, dbeT = sp.symbols(f'dal{s_} dbe{s_} dal{s_}T dbe{s_}T'); unk += [dal, dbe] + ([dalT, dbeT] if hasT else [])
        add(sp.expand(expr[s_] - (al0 * I1 + (dal + (TT * dalT if hasT else 0)) * I0 + dbe + (TT * dbeT if hasT else 0))))
    M, rhs = sp.linear_eq_to_matrix(eqs, unk); rk, rka = M.rank(), sp.Matrix.hstack(M, rhs).rank()
    return rk, rka, len(unk), len(eqs), (sp.linsolve((M, rhs), unk), unk) if rk == rka else None

if __name__ == "__main__":
    t0 = time.time()
    # P1: independent recomputation of four columns at spins 1-5 with code/wkb.py
    V0 = sp.Rational(N0, 2) * X * u/(1+u) + (sp.Rational(N0 + 2, 2) * Y + sp.Rational(1, 4)) * u/(1+u)**2
    Sm = wronskian_S(V0, 5)
    mine = {'m': {s_: sp.expand(sp.diff(Sm[s_], msym).subs(msym, N0)) for s_ in (1, 3, 5)}}
    c = sp.Symbol('c')
    for name, f in (('b1', u/(1+u)), ('b2', u/(1+u)**2), ('b3', u/(1+u)**3)):
        S = wronskian_S(V0 + c * f, 5, msym)
        mine[name] = {s_: sp.expand(sp.diff(S[s_], c).subs({c: 0, msym: N0})) for s_ in (1, 3, 5)}
    cols = {nm: load_col(f"n-1_{nm}") for nm in ('m', 'b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'a1', 'a2', 'L1', 'L2', 'Z1', 'Z2')}
    ok1 = all(sp.expand(mine[nm][s_] - cols[nm][s_]) == 0 for nm in mine for s_ in (1, 3, 5))
    require(ok1, "P1: pinned columns of delta m, u/(1+u), u/(1+u)^2, u/(1+u)^3 at spins 1, 3, 5 equal an independent recomputation with code/wkb.py")
    # P2: base
    base = load_col("base_S")
    ok2 = True
    for s_ in (3, 5, 7):
        pc = sp.expand(ALPHA[s_].subs(n, N0) * IPC[s_].subs(n, N0) + BETA[s_])
        ok2 &= sp.expand(base[s_] - pc) == 0
    require(ok2, "P2: pinned base S_3, S_5, S_7 equal alpha_s I_s^pc + beta_s of the paper at n = -1 (with the Bernoulli beta_s)")
    for sol in (3, 1):
        Ic = cyl(sol, 10); I0 = sp.expand(sum(sp.limit(c_, t, sp.oo) * X**a * Y**b for (a, b), c_ in sp.Poly(Ic, X, Y).terms()))
        Pm = sp.Poly(base[9] - base[9].subs({X: 0, Y: 0}), X, Y); Pi = sp.Poly(I0 - I0.subs({X: 0, Y: 0}), X, Y)
        ratios = {mon: sp.nsimplify(cS / Pi.coeff_monomial(X**mon[0] * Y**mon[1])) for mon, cS in Pm.terms()}
        al9 = sp.Rational(json.load(open(os.path.join(FO, f"alpha9_sol{sol}.json")))['alpha9_0'])
        require(len(set(ratios.values())) == 1 and list(ratios.values())[0] == al9 and al9 == sp.Rational(-1, 1440),
                f"P2: spin-9 correspondence, Solution {sol}: the base S_9 is alpha_9 x the weight-10 table at t = infinity over all 14 momentum monomials, alpha_9 = -1/1440")
    # P3: positive control
    ctrl = {}
    for s_ in (3, 5, 7):
        I0 = sp.expand(LAM[s_] * IPC[s_].subs(n, N0)); I1 = sp.expand(LAM[s_] * sp.diff(IPC[s_], n).subs(n, N0))
        ctrl[s_] = (ALPHA[s_].subs(n, N0) / LAM[s_], I0, I1)
    rk, rka, nu, ne, res = solve(cols, ctrl, ['b1', 'b2'], (3, 5, 7), spin1=False)
    expect = {'dm': 1, 'c1_b1': sp.Rational(1, 2), 'c2_b2': sp.Rational(1, 2)}
    got = {}
    if res:
        for s__ in res[0]:
            got = {str(u_): v_ for u_, v_ in zip(res[1], s__) if v_ != 0}
    require(res is not None and got == expect, f"P3: POSITIVE CONTROL -- the paperclip's own n-deformation is recovered exactly (delta m = 1, delta A = X/2, delta B = Y/2, nothing else): got {got}")
    # P4: the classes, spins 1-9
    for sol in (3, 1):
        data = first_order(sol, (3, 5, 7, 9))
        for label, terms in (("the two LVZ terms", ['b1', 'b2']),
                             ("the (1+u)-library with growth terms", ['b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'a1', 'a2']),
                             ("the library with the moving-exponent log columns", ['b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'a1', 'a2', 'L1', 'L2', 'Z1', 'Z2'])):
            rk, rka, nu, ne, res = solve(cols, data, terms, (3, 5, 7, 9))
            require(res is None and rka == rk + 1, f"P4: Solution {sol}, spins 1-9: {label} is INCONSISTENT ({ne} equations, {nu} unknowns, rank {rk}, augmented {rka})")
        data7 = first_order(sol, (3, 5, 7))
        for v0 in ("3_1", "9_1", "1_2", "2_1"):
            fp = dict(cols); fp[f'D{v0}'] = load_col(f"fp_v{v0}_D" if v0 != "2_1" else "fp_v2_1_D5"); fp[f'E{v0}'] = load_col(f"fp_v{v0}_E" if v0 != "2_1" else "fp_v2_1_E5")
            spins = (3, 5, 7) if v0 != "2_1" else (3, 5)
            rk, rka, nu, ne, res = solve(fp, data7, ['b1', 'b2', 'b3', 'b4', f'D{v0}', f'E{v0}'], spins)
            require(res is None and rka == rk + 1, f"P4: Solution {sol}: a fixed fourth point at u = -{v0.replace('_', '/')} (spins {spins}) is INCONSISTENT ({ne} equations, {nu} unknowns, rank {rk}, augmented {rka})")
    # P5: negative controls
    ctrl_bad = dict(ctrl); s_ = 5; al, I0, I1 = ctrl[s_]; ctrl_bad[s_] = (al, I0, sp.expand(I1 + X))
    rk, rka, nu, ne, res = solve(cols, ctrl_bad, ['b1', 'b2'], (3, 5, 7), spin1=False)
    require(res is None, "P5: negative control -- a tampered data coordinate makes the positive control inconsistent (guard live)")
    rk, rka, nu, ne, res = solve(cols, ctrl, ['b1', 'b2'], (3, 5, 7), spin1=False)
    got = {str(u_): v_ for s__ in res[0] for u_, v_ in zip(res[1], s__) if v_ != 0} if res else {}
    require(got != {'dm': 2, 'c1_b1': sp.Rational(1, 2), 'c2_b2': sp.Rational(1, 2)}, "P5: negative control -- the positive control does not accept delta m = 2 (guard live)")
    print(f"\n{time.time()-t0:.0f}s")
    if fails:
        print(f"FIRST ORDER: {len(fails)} FAILURE(S)"); sys.exit(1)
    print("FIRST ORDER CERTIFIED (pinned response columns cross-checked with code/wkb.py; the paperclip correspondence holds at "
          "spin 9 with alpha_9 = -1/1440; the pipeline recovers the paperclip's own deformation exactly; with spins 1-9 the two-term, "
          "(1+u)-library and moving-exponent classes and, with spins 1-7, a fixed fourth point at u = -3, -9, -1/2, -2 are all "
          "inconsistent with the first-order deformation of Solutions 3 and 1)")
