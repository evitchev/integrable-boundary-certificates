"""THE FIRST-ORDER DEFORMATION OF THE c = -2 ANCHOR IS BILINEAR IN THE TWO SECTOR DETERMINANTS (Solution 3 at t = 3).
Exact data: the top-normalized vacuum charges e^_(2k-1)(t) of Solution 3 (rational in t; code/cyl_first_order tables through W = 10 and
the hold-out-verified vev profiles at W = 12, 14), their derivative at t = 3.  Sector determinants D_S(E) = Gamma(1+beta_S)/Gamma(a_S),
a_S = (1+beta_S)/2 - E/4 (Masoero-Ruzza (5.32)), beta_X^2 = -4X, beta_Y^2 = 1 - 4Y; the vacuum e^_(2k-1) = (-2)^k [F_k(D_X) + F_k(D_Y)] is
the E^(1-2k) coefficient of log(D_X D_Y) in units s_k = 2^(5k-3)/(k(2k-1)) (both re-verified here).
RESULTS (all exact):
 (1) TOP LAYER: the degree-k part of d/dt e^_(2k-1)|_(t=3) is  sum_{i=1}^{k-1} k(2k-1)/(4 i (k-i)) X^(k-i) Y^i  -- read off at k = 2..5 and
     PREDICTED at k = 6, 7 (spins 11, 13): 9 hold-out coefficients, all confirmed.  Generating function: (E/32) log(1+16X/E^2) log(1+16Y/E^2)
     = 2E (dL_X/dE)(dL_Y/dE), L_S the classical top layer of log D_S: a PRODUCT of the two sectors' classical level densities.
 (2) ALL MIXED COEFFICIENTS: with the exact quantum level densities B_S(E) := [d/dE log D_S]_(even in beta_S) (the beta-odd part carries the
     sector's own even-E coefficients and cannot enter a polynomial), the E^(1-2k) coefficients of  2E B_X(E) B_Y(E), converted to charge
     units and normalized on X^k, reproduce EVERY mixed coefficient X^a Y^b (a, b >= 1) of d/dt e^ at every degree, spins 3..13 -- 56 mixed
     coefficients of which only the top layer at k <= 5 (10 numbers) entered the derivation: 46 parameter-free predictions confirmed.
     REFUTED AS WRITTEN (not an equivalent form -- Codex, reviews of 30738b9 and queued batch 34; sec. 38(am) correction 1):
     D(t; E) = D_X D_Y +
     2 (t-3) E D_X'(E) D_Y'(E) + ... with the root law delta E_n^X = -2(t-3) E_n^X sum_k 1/(E_n^X - E_k^Y).  Passing to an analytic
     determinant discards the beta-even PROJECTION that makes the identity above work: unprojected, the E^-3 coefficient of
     2E R_X R_Y differs by (2/3) beta_X beta_Y (beta_X^2 + beta_Y^2 - 2), a mixed non-polynomial term no single-sector term cancels.
     So the displayed determinant form is REFUTED as written, not merely unproven; what this script certifies is the projected
     asymptotic identity, and an analytic determinant realizing it is an open construction.
 (3) SINGLE-SECTOR REMAINDER (pure-X plus pure-Y polynomials): Y-sector top = (3/16) d/dY e^_Y (a momentum shift of the Y sector, none for
     X); both sectors then carry the SAME weight-2 term with top -(1/960) E d/dE d_S^2 [log D_S] = (1/120) E (d_S B_S)^2 (identical at top,
     different below); the remainder's next layers are again common to the two sectors (quadratic-in-k tops); its closed form is OPEN.
     Pinned negative: the 11-operator basis {d_S, EdE d_S, E^-2, E^-1 B, d_S^2, EdE d_S^2, E(d_S B)^2, E^-1 (d_S A)^2, (d_S A)(d_S B),
     E^-2 d_S, E^-1 d_S B} on the exact sector series (22 unknowns) does NOT close (Groebner [1] through spin 11 or hold-out failure).
 CONTROLS (must FIRE): (c1) coefficient 1 instead of 2 in the bilinear; (c2) E^3 instead of E; (c3) the product with the full level densities
     (beta-odd parts kept) is not a polynomial in X, Y.
Usage: t3_bilinear.py   (exit 0 iff every check passes and every control fires)"""
import sys, os, json
ROOT = os.environ.get('IB_ROOT', os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, 'code')); os.chdir(os.path.join(ROOT, 'code'))
import sympy as sp
from cyl_first_order import cyl, t, X, Y
KM = 7; K = 2 * KM + 2
E, a, eps = sp.symbols('E a eps'); bX, bY = sp.symbols('beta_X beta_Y')
failures = []
def check(name, ok, detail=''):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  {detail}" if detail else ''), flush=True)
    if not ok: failures.append(name)
def sol3_charge(W):
    if W <= 10: return sp.expand(cyl(3, W))
    f = {12: 'results/lab/anchor11/vev_profile_sol3_w12.json', 14: 'results/lab/anchor13/vev_profile_sol3_w14.json'}[W]
    d = json.load(open(os.path.join(ROOT, f))); assert d['solution'] == 3 and d['weight'] == W and len(d['hold_outs']) >= 4
    return sp.expand(sum(sp.sympify(v) * X**int(kk.split(',')[0]) * Y**int(kk.split(',')[1]) for kk, v in d['vev'].items()))
def dlogD_series(b):
    """psi(a)/4 = d/dE log D as a series in 1/E (the constant log(-E/4) of log a dropped)"""
    psi = sp.log(1 - 2 * (1 + b) * eps) - 1 / (2 * a) - sum(sp.bernoulli(2 * k) / (2 * k * a**(2 * k)) for k in range(1, K // 2 + 3))
    ser = sp.series((psi / 4).subs(a, (1 + b) / 2 - 1 / (4 * eps)), eps, 0, K + 3).removeO()
    return sp.expand(sp.expand(ser).subs(eps, 1 / E))
def logD_series(b):
    st = -((a - sp.Rational(1, 2)) * sp.log(1 - 2 * (1 + b) * eps) - a) - sum(sp.bernoulli(2 * k) / (2 * k * (2 * k - 1) * a**(2 * k - 1)) for k in range(1, KM + 3))
    ser = sp.series(st.subs(a, (1 + b) / 2 - 1 / (4 * eps)), eps, 0, K + 3).removeO()
    return sp.expand(sp.expand(ser).subs(eps, 1 / E))
ev = lambda ex, b: sp.expand((ex + ex.subs(b, -b)) / 2)
sk = {k: sp.Rational(2**(5 * k - 3), k * (2 * k - 1)) for k in range(1, KM + 1)}
D = sp.Symbol('Delta'); z = (1 + sp.sqrt(1 + 8 * D)) / 2
Fk = {k: sp.expand(sp.simplify(sp.bernoulli(2 * k, z) / 2**k)) for k in range(1, KM + 1)}
lX, lY = logD_series(bX), logD_series(bY); dX, dY = dlogD_series(bX), dlogD_series(bY)
check("s_k = 2^(5k-3)/(k(2k-1)) is the vacuum normalization (E^(1-2k) coefficient of log D_S = s_k F_k(Delta_S)), k = 1..7",
      all(sp.simplify(lX.coeff(E, -(2 * k - 1)) / sp.expand(Fk[k].subs(D, (bX**2 - 1) / 8)) - sk[k]) == 0 for k in range(1, KM + 1)))
check("parity: [log D]_even-in-beta has odd powers of E only; [d/dE log D]_even has even powers only",
      all(ev(lX, bX).coeff(E, -j) == 0 for j in range(2, K, 2)) and all(ev(dX, bX).coeff(E, -j) == 0 for j in range(1, K, 2)))
def to_XY(expr):
    p = sp.Poly(sp.expand(expr), bX, bY); out = 0
    for (i, j), c in p.terms():
        if i % 2 or j % 2: return None
        out += c * (-4 * X)**(i // 2) * (1 - 4 * Y)**(j // 2)
    return sp.expand(out)
targets = {}
for k in range(2, KM + 1):
    e = sol3_charge(2 * k); p = sp.Poly(e, X, Y, domain=sp.QQ.frac_field(t))
    top = sp.QQ.frac_field(t).to_sympy(p.coeff_monomial(X**k)); en = sp.cancel(e / top)
    e0 = sp.expand(en.subs(t, 3)); d1 = sp.expand(sp.diff(en, t).subs(t, 3)); targets[k] = (e0, d1)
check("vacuum e^_(2k-1)(t=3) = (-2)^k [F_k(-X/2-1/8) + F_k(-Y/2)], spins 3-13",
      all(sp.expand(targets[k][0] - (-2)**k * (Fk[k].subs(D, -X / 2 - sp.Rational(1, 8)) + Fk[k].subs(D, -Y / 2))) == 0 for k in range(2, KM + 1)))
# (1) top layer
tau = {k: {m_: c_ for m_, c_ in sp.Poly(targets[k][1], X, Y).terms() if sum(m_) == k} for k in range(2, KM + 1)}
closed = lambda k, i: sp.Rational(k * (2 * k - 1), 4 * i * (k - i))
check("(1) top layer tau_{k,i} = k(2k-1)/(4 i (k-i)) at k = 2..5 (spins 3-9, the fit range)",
      all(tau[k].get((k - i, i), 0) == closed(k, i) for k in range(2, 6) for i in range(1, k)) and all(tau[k].get((k, 0), 0) == 0 and tau[k].get((0, k), 0) == 0 for k in range(2, 6)))
check("(1) HOLD-OUT: the same closed form at k = 6, 7 (spins 11, 13; 9 predicted coefficients)",
      all(tau[k].get((k - i, i), 0) == closed(k, i) for k in (6, 7) for i in range(1, k)) and all(tau[k].get((k, 0), 0) == 0 and tau[k].get((0, k), 0) == 0 for k in (6, 7)))
# generating function of the top layer: (E/32) log(1+16X/E^2) log(1+16Y/E^2) -> E^(1-2k) coefficient x (-2)^k / s_k
u = sp.Symbol('u'); G = sp.expand(sp.series(sp.log(1 + 16 * X * u) * sp.log(1 + 16 * Y * u) / 32, u, 0, KM + 1).removeO())   # u = E^-2; E x u^k -> E^(1-2k)
check("(1) generating function (E/32) log(1+16X/E^2) log(1+16Y/E^2) reproduces the top layer, k = 2..7",
      all(sp.expand((-2)**k / sk[k] * G.coeff(u, k) - sum(closed(k, i) * X**(k - i) * Y**i for i in range(1, k))) == 0 for k in range(2, KM + 1)))
# (2) all mixed coefficients from the bilinear 2E B_X B_Y
def mixed_test(expr, label, expect_pass):
    ok = True; nmixed = 0; detail = []
    for k in range(2, KM + 1):
        c = to_XY(ev(ev(expr.coeff(E, -(2 * k - 1)), bX), bY))
        if c is None: ok = False; detail.append(f"spin {2*k-1}: not polynomial"); continue
        M = sp.expand((-2)**k / sk[k] * c); e0, d1 = targets[k]
        pred = sp.expand(M - e0 * sp.Poly(M, X, Y).coeff_monomial(X**k))
        diff = sp.Poly(sp.expand(d1 - pred), X, Y)
        mixed = [(m_, c_) for m_, c_ in diff.terms() if m_[0] > 0 and m_[1] > 0]
        nm = len([m_ for m_, c_ in sp.Poly(d1, X, Y).terms() if m_[0] > 0 and m_[1] > 0]); nmixed += nm
        if mixed: ok = False; detail.append(f"spin {2*k-1}: {len(mixed)} mixed residuals")
    check(f"(2) {label}: mixed part of d/dt e^ reproduced at spins 3-13 ({nmixed} mixed coefficients)" if expect_pass else f"CONTROL {label} fires",
          ok if expect_pass else not ok, '; '.join(detail))
    return ok
dXe, dYe = ev(dX, bX), ev(dY, bY)
mixed_test(sp.expand(2 * E * dXe * dYe), "bilinear 2E B_X B_Y (exact quantum level densities, beta-even parts)", True)
mixed_test(sp.expand(E * dXe * dYe), "(c1) coefficient 1", False)
mixed_test(sp.expand(2 * E**3 * dXe * dYe), "(c2) E^3 instead of E", False)
full = sp.expand(2 * E * dX * dY)
check("CONTROL (c3) the product with the full level densities (beta-odd parts kept) is NOT polynomial in X, Y",
      any(to_XY(full.coeff(E, -(2 * k - 1))) is None for k in range(2, KM + 1)))
# (3) single-sector remainder: tops
def series_poly(expr_E, b, bsub, powers):
    out = {}
    for j in powers:
        c = ev(expr_E.coeff(E, -j), b); p = sp.Poly(c, b); out[j] = sp.expand(sum(cc * bsub**(ii // 2) for (ii,), cc in p.terms()))
    return out
A = {'X': sum(v * E**(-j) for j, v in series_poly(lX, bX, -4 * X, range(1, K, 2)).items()), 'Y': sum(v * E**(-j) for j, v in series_poly(lY, bY, 1 - 4 * Y, range(1, K, 2)).items())}
B = {'X': sum(v * E**(-j) for j, v in series_poly(dX, bX, -4 * X, range(2, K, 2)).items()), 'Y': sum(v * E**(-j) for j, v in series_poly(dY, bY, 1 - 4 * Y, range(2, K, 2)).items())}
def cu(expr, k): return sp.expand((-2)**k / sk[k] * sp.expand(expr).coeff(E, -(2 * k - 1)))
resid = {}
bil = sp.expand(2 * E * dXe * dYe)
for k in range(2, KM + 1):
    M = sp.expand((-2)**k / sk[k] * to_XY(bil.coeff(E, -(2 * k - 1)))); e0, d1 = targets[k]
    resid[k] = sp.expand(d1 - sp.expand(M - e0 * sp.Poly(M, X, Y).coeff_monomial(X**k)))
okY = all(sp.Poly(resid[k], X, Y).coeff_monomial(Y**(k - 1)) == sp.Poly(cu(sp.Rational(3, 16) * sp.diff(A['Y'], Y), k), X, Y).coeff_monomial(Y**(k - 1)) and sp.Poly(resid[k], X, Y).coeff_monomial(X**(k - 1)) == 0 for k in range(2, KM + 1))
check("(3) single-sector remainder: Y top = (3/16) d/dY e^_Y at every spin, X has no degree-(k-1) term", okY)
w2 = {}
for k in range(3, KM + 1):
    ra = cu(-E * sp.diff(sp.diff(A['X'], X, 2), E) / 960, k); rb = cu(E * sp.diff(B['X'], X)**2 / 120, k)
    w2[k] = (sp.Poly(resid[k], X, Y).coeff_monomial(X**(k - 2)), sp.Poly(resid[k], X, Y).coeff_monomial(Y**(k - 2)) - cu(sp.Rational(3, 16) * sp.diff(A['Y'], Y), k).coeff(Y**(k - 2)) if False else sp.Poly(sp.expand(resid[k] - cu(sp.Rational(3, 16) * sp.diff(A['Y'], Y), k)), X, Y).coeff_monomial(Y**(k - 2)),
             sp.Poly(ra, X, Y).coeff_monomial(X**(k - 2)), sp.Poly(rb, X, Y).coeff_monomial(X**(k - 2)))
check("(3) weight-2 top: X-part and (Y-part minus the shift) both = k(2k-1)(k-1)/960 S^(k-2) = -(1/960) E d/dE d_S^2 log D_S = (1/120) E (d_S B_S)^2 at top, k = 3..7",
      all(w2[k][0] == w2[k][1] == w2[k][2] == w2[k][3] == sp.Rational(k * (2 * k - 1) * (k - 1), 960) for k in range(3, KM + 1)))
print("remainder after the Y shift and the weight-2 term (form (a): -(1/960) E d_E d_S^2 log D_S in both sectors):")
for k in range(4, KM + 1):
    r = sp.expand(resid[k] - cu(sp.Rational(3, 16) * sp.diff(A['Y'], Y), k) - cu(-E * sp.diff(sp.diff(A['X'], X, 2) + sp.diff(A['Y'], Y, 2), E) / 960, k))
    r = sp.expand(r - targets[k][0] * sp.Poly(r, X, Y).coeff_monomial(X**k)); p = sp.Poly(r, X, Y)
    print(f"   spin {2*k-1}: X layers {dict(sorted({m_[0]: c_ for m_, c_ in p.terms() if m_[1] == 0 and m_[0] > 0}.items(), reverse=True))}; Y layers {dict(sorted({m_[1]: c_ for m_, c_ in p.terms() if m_[0] == 0 and m_[1] > 0}.items(), reverse=True))}")
print("SUMMARY:", "ALL CHECKS PASS, ALL CONTROLS FIRE" if not failures else f"FAILURES: {failures}")
sys.exit(0 if not failures else 1)
