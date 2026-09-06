"""Heun-type fourth singular point at m = -1: rationality (no-transcendental) conditions of the
WKB coefficients S_1, S_3, S_5, S_7, reconstructed as COMPLETE rational functions of the pole
position v (prefactors included) from exact evaluations at many rational v, with hold-out
verification, FAIL-CLOSED (Codex round 21).  Writes results/lab/heun/heun_rationality.json.

Potential: V = A u/(1+u) + B u/(1+u)^2 + D u/(u+v) + E u/(u+v)^2, leading symbol kappa^2/(1+u).
At each v the T-coefficient of S_k is an exact polynomial in the coefficient symbols; every
polynomial coefficient is reconstructed in v with a degree bound DEG (default 8) that must be
strictly respected by the fit (a coefficient needing degree > DEG is a failure), and verified
on HOLD hold-out positions never used in the fit.  Classification: exact finite-sample
evidence with verified low-degree reconstructions; the degree bound is an input, not derived.

Stages (each imposes the relations found at the previous stage):
  spin 1: full T-coefficient;   spin 3 (E = 2D(v-1));   spin 5 (+ spin-3 relation, B eliminated);
  spin 7 (+ spin-5 linear relation, momentum-dependent branch).  Runtime ~1 h (spin 7 dominates).
"""
import sys, json, time, os
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(os.path.dirname(HERE), "code"))
import sympy as sp
from fractions import Fraction as F
import cyl_wkb_heun as H
from cyl_wkb_heun import u, v, T, A, B, D, E
from invariant_engine import nullspace_cols

DEG = 8          # fixed degree bound (recorded in the output JSON)
FIT = [F(a, b) for a, b in [(3,1),(7,3),(5,2),(9,1),(4,1),(6,1),(11,4),(13,5),(17,3),(8,1),(12,5),(7,2),(19,7),(10,1),(23,9),(15,2),(11,3),(5,1),(29,13),(14,3)]]
HOLD = [F(a, b) for a, b in [(21,11),(16,3),(31,7),(2,1)]]
FIT7 = FIT; HOLD7 = HOLD   # spin 7 on ALL positions: the fitter needs 2(deg+1) samples, and 8 positions capped the degree at 3 (run of 2026-09-03 15:09: D^2 and D^4 unfittable)
vs = sp.Symbol('v')
failures = []

def tcoeff(rk, k, vq):
    rk = sp.together(rk)
    num, den = sp.fraction(rk)
    dpoly = sp.Poly(den, u); b = d = 0
    while True:
        q_, rem = sp.div(dpoly, sp.Poly(1 + u, u))
        if rem.is_zero and dpoly.degree() > 0: dpoly = q_; b += 1
        else: break
    while True:
        q_, rem = sp.div(dpoly, sp.Poly(u + vq, u))
        if rem.is_zero and dpoly.degree() > 0: dpoly = q_; d += 1
        else: break
    if dpoly.degree() != 0:
        raise ValueError("denominator")
    const = dpoly.coeffs()[0]
    npoly = sp.Poly(sp.expand(num), u)
    total = sum(cf * H.euler(a, b - sp.Rational(k, 2), d).subs(v, vq) for (a,), cf in zip(npoly.monoms(), npoly.coeffs())) / const
    pl = sp.Poly(sp.expand(total), T)
    return sp.expand(pl.coeff_monomial(T)) if pl.degree() >= 1 else sp.Integer(0)

def fit_rational(samples, deg):
    ts = sorted(samples)
    cols = [[t ** k for t in ts] for k in range(deg + 1)] + [[-samples[t] * t ** k for t in ts] for k in range(deg + 1)]
    ker = nullspace_cols(cols, len(ts))
    return None if not ker else (list(ker[0][:deg + 1]), list(ker[0][deg + 1:]))

def ev(c, t): return sum(ci * t ** k for k, ci in enumerate(c))

def reconstruct(name, samples, hold):
    """rational function of v of degree <= DEG fitting all samples and verified on all hold-outs; fail-closed"""
    for deg in range(0, DEG + 1):
        if len(samples) < 2 * (deg + 1):
            break
        r = fit_rational(samples, deg)
        if r is None: continue
        P_, Q_ = r
        if all(ev(Q_, t) != 0 and ev(P_, t) / ev(Q_, t) == samples[t] for t in samples) and all(ev(Q_, t) != 0 and ev(P_, t) / ev(Q_, t) == hold[t] for t in hold):
            num = sum(sp.Rational(c.numerator, c.denominator) * vs ** k for k, c in enumerate(P_)); den = sum(sp.Rational(c.numerator, c.denominator) * vs ** k for k, c in enumerate(Q_))
            expr = sp.factor(sp.cancel(num / den))
            print(f"PASS  {name} = {expr}   [deg {deg} <= {DEG}, {len(samples)} fit, {len(hold)} hold-out]", flush=True)
            return expr
    print(f"FAIL  {name}: no rational function of degree <= min({DEG}, {len(samples)//2 - 1}) (degree bound, sample count) fits the {len(samples)} samples and verifies on {len(hold)} hold-outs", flush=True)
    failures.append(name); return None

def poly_coeffs(expr, gens):
    pl = sp.Poly(sp.expand(expr), *gens)
    return {mon: F(str(c)) for mon, c in zip(pl.monoms(), pl.coeffs())}

out = {"deg_bound": DEG, "fit_positions": [str(t) for t in FIT], "holdout_positions": [str(t) for t in HOLD], "stages": {}}
t0 = time.time()
# ---- spin 1 (full), then spin 3 with E = 2D(v-1)
V = A*u/(1+u) + B*u/(1+u)**2 + D*u/(u+v) + E*u/(u+v)**2
r = H.riccati(V, 3)
for k, gens, tag in ((1, (A, B, D, E), "spin1"), (3, (A, B, D, E), "spin3_full")):
    samples = {}
    for vv in FIT + HOLD:
        vq = sp.Rational(vv.numerator, vv.denominator)
        samples[vv] = poly_coeffs(tcoeff(r[k].subs(v, vq), k, vq), gens)
    mons = sorted({m_ for s_ in samples.values() for m_ in s_})
    rec = {}
    for m_ in mons:
        name = f"{tag} T-coefficient of " + "*".join(f"{g}^{e}" for g, e in zip(gens, m_) if e) if any(m_) else f"{tag} T-coefficient constant"
        rec[str(m_)] = str(reconstruct(name, {t: samples[t].get(m_, F(0)) for t in FIT}, {t: samples[t].get(m_, F(0)) for t in HOLD}))
    out["stages"][tag] = {"gens": [str(g) for g in gens], "coefficients": rec}
print(f"spin 1 and spin 3 (full) done ({time.time()-t0:.0f}s)", flush=True)
# ---- spin 3 with E = 2D(v-1): coefficients of A, B, D, 1 after dividing by D
V3 = A*u/(1+u) + B*u/(1+u)**2 + D*u/(u+v) + 2*D*(v-1)*u/(u+v)**2
r3 = H.riccati(V3, 3)
samples = {}
for vv in FIT + HOLD:
    vq = sp.Rational(vv.numerator, vv.denominator)
    tc = tcoeff(r3[3].subs(v, vq), 3, vq)
    q, rem = sp.div(sp.Poly(tc, D), sp.Poly(D, D))
    if not rem.is_zero:
        print(f"FAIL  spin-3 T-coefficient not divisible by D at v = {vv}", flush=True); failures.append("spin3 divisibility"); break
    samples[vv] = poly_coeffs(q.as_expr(), (A, B, D))
if "spin3 divisibility" not in failures:
    rec = {}
    for m_, nm in (((1,0,0), "A"), ((0,1,0), "B"), ((0,0,1), "D"), ((0,0,0), "1")):
        rec[nm] = reconstruct(f"spin3/D coefficient of {nm}", {t: samples[t].get(m_, F(0)) for t in FIT}, {t: samples[t].get(m_, F(0)) for t in HOLD})
    out["stages"]["spin3_reduced"] = {k_: str(v_) for k_, v_ in rec.items()}
    cA, cB, cD, c0 = rec["A"], rec["B"], rec["D"], rec["1"]
print(f"spin 3 reduced done ({time.time()-t0:.0f}s)", flush=True)
# ---- spin 5 with B eliminated: full T-coefficient as a polynomial in A, D (coefficients reconstructed, prefactor included)
if all(x is not None for x in (cA, cB, cD, c0)):
    samples = {}
    for vv in FIT + HOLD:
        vq = sp.Rational(vv.numerator, vv.denominator)
        Bsol = -(cA*A + cD*D + c0).subs(vs, vq) / cB.subs(vs, vq)
        V5 = A*u/(1+u) + Bsol*u/(1+u)**2 + D*u/(u+vq) + 2*D*(vq-1)*u/(u+vq)**2
        r5 = H.riccati(V5, 5)
        samples[vv] = poly_coeffs(tcoeff(r5[5], 5, vq), (A, D))
    mons = sorted({m_ for s_ in samples.values() for m_ in s_})
    rec = {}
    for m_ in mons:
        rec[str(m_)] = str(reconstruct(f"spin5 T-coefficient of A^{m_[0]} D^{m_[1]}", {t: samples[t].get(m_, F(0)) for t in FIT}, {t: samples[t].get(m_, F(0)) for t in HOLD}))
    out["stages"]["spin5_full"] = {"gens": ["A", "D"], "coefficients": rec}
    # factorization of the reconstructed polynomial in (A, D) over Q(v)
    poly5 = sum(sp.sympify(c) * A**m_[0] * D**m_[1] for m_, c in ((eval(k_), c_) for k_, c_ in rec.items()) if c != "None")
    fac5 = sp.factor(sp.together(poly5))
    out["stages"]["spin5_factored"] = str(fac5)
    out["stages"]["spin5_prefactor_zeros"] = [str(z_) for z_ in sp.solve(sp.numer(sp.together(sp.Poly(sp.numer(sp.together(poly5)), A, D).LC())), vs)]
    print(f"PASS  spin 5 reconstructed T-coefficient factors as {fac5}", flush=True)
    print(f"spin 5 done ({time.time()-t0:.0f}s)", flush=True)
    # ---- spin 7 on the momentum-dependent branch (linear spin-5 factor solved for A), full polynomial in D
    lin = [f_ for f_, e_ in sp.factor_list(sp.numer(sp.together(poly5)))[1] if f_.has(A)]
    if len(lin) == 1:
        Asol = sp.solve(lin[0], A)[0]
        samples = {}
        for vv in FIT7 + HOLD7:
            vq = sp.Rational(vv.numerator, vv.denominator)
            Bsol = -(cA*Asol.subs(vs, vq) + cD*D + c0).subs(vs, vq) / cB.subs(vs, vq)
            V7 = Asol.subs(vs, vq)*u/(1+u) + Bsol*u/(1+u)**2 + D*u/(u+vq) + 2*D*(vq-1)*u/(u+vq)**2
            r7 = H.riccati(V7, 7)
            samples[vv] = poly_coeffs(tcoeff(r7[7], 7, vq), (D,))
            print(f"      spin 7 at v = {vv} ({time.time()-t0:.0f}s)", flush=True)
        mons = sorted({m_ for s_ in samples.values() for m_ in s_})
        rec = {}
        for m_ in mons:
            rec[str(m_)] = str(reconstruct(f"spin7 T-coefficient of D^{m_[0]}", {t: samples[t].get(m_, F(0)) for t in FIT7}, {t: samples[t].get(m_, F(0)) for t in HOLD7}))
        out["stages"]["spin7_full"] = {"gens": ["D"], "coefficients": rec}
        poly7 = sum(sp.sympify(c) * D**m_[0] for m_, c in ((eval(k_), c_) for k_, c_ in rec.items()) if c != "None")
        fac7 = sp.factor(sp.together(poly7))
        out["stages"]["spin7_factored"] = str(fac7)
        lc7 = sp.Poly(sp.numer(sp.together(poly7)), D).LC()
        pref_zeros = sp.solve(sp.numer(sp.together(lc7)), vs)
        out["stages"]["spin7_prefactor_zeros"] = [str(z_) for z_ in pref_zeros]
        print(f"PASS  spin 7 reconstructed T-coefficient factors as {fac7}; leading-coefficient zeros in v: {pref_zeros}", flush=True)
    else:
        print("FAIL  spin-5 factorization has no unique A-linear factor", flush=True); failures.append("spin5 factor shape")
os.makedirs(os.path.join(os.path.dirname(HERE), "results", "lab", "heun"), exist_ok=True)
out["failures"] = failures
json.dump(out, open(os.path.join(os.path.dirname(HERE), "results", "lab", "heun", "heun_rationality.json"), "w"), indent=1)
print(f"\n{time.time()-t0:.0f}s")
if failures:
    print(f"HEUN RATIONALITY: {len(failures)} FAILURE(S)"); sys.exit(1)
print("HEUN RATIONALITY RECONSTRUCTED (exact finite-sample evidence: rational expressions incl. prefactors fitted under an ASSERTED degree bound DEG and verified on hold-outs; not a proof for all v)")
