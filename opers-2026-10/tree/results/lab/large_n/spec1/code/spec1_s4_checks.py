"""SPEC1 S4: machine checks of the analytic steps (alpha)-(gamma) in the general-family proof (SPEC1_REPORT.md sec. S4).  NEW FILE.
Usage: spec1_s4_checks.py S4
General family: c = sigma (a+1), n = sigma a, D = d + a, q_G = c/D, q_O = c/(d-1), p = q_G/sigma = (a+1)/(d+a), d >= 2.
 (C1) exact: p < 1 iff d > 1; the E-term corrections to theta S are x^(n - (d-1) q_O) = x^(-sigma) (ODE) and
      x^(n - (D-1) q_G) = x^(q_G - sigma) = x^(sigma (p-1)) (Gamma form): negative powers, so E enters neither the exponent nor beta.
 (C2) exact: Laplace scaling.  tau* = x^(q_O)/sigma solves the saddle equation p x^(q_G) sigma^p tau^(p-1)/q_G = 1; with tau = tau* v the
      phase tau - X^(q_G)/q_G equals tau* h(v), h(v) = v - v^p/p, and tau* h(1) = -x^(q_O)/q_O; h''(1) = 1 - p.
 (C3) the descent inequality on the vertical line v = 1 + i y:  Re h(1 + i y) - h(1) = -R(y)/p,
      R(y) = cos(p t)/cos(t)^p - 1, t = arctan y; R > 0 for y != 0 (strict concavity of log cos), R = p(1-p) y^2/2 + O(y^4),
      R >= cos(p pi/2) |y|^p - 1.  Exact series + identity; dps-30 grid of R/min(y^2, |y|^p) at the six p of the S3 points (info).
 (C4) Laplace main term on the vertical line: sqrt(lam)/(2 pi) int e^(lam (h(1+iy) - h(1))) (1+iy)^kappa dy -> (2 pi (1-p))^(-1/2),
      error O(1/lam): quadrature at lam = 50, 200, 800 for p = 4/5, kappa = 0.3 - 0.2 i (info: rate check)."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import spec1_lib as SL  # noqa: E402
from spec1_lib import sp, mp  # noqa: E402

FAIL = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""), flush=True)
    if not ok:
        FAIL.append(name)


def c1_c2():
    a, d, sig, x, v = sp.symbols("a d sigma x v", positive=True)
    c = sig * (a + 1)
    n = sig * a
    D = d + a
    qG, qO = c / D, c / (d - 1)
    p = sp.simplify(qG / sig)
    check("(C1) p = q_G/sigma = (a+1)/(d+a); 1 - p = (d-1)/(d+a)", sp.simplify(p - (a + 1) / (d + a)) == 0
          and sp.simplify(1 - p - (d - 1) / (d + a)) == 0)
    check("(C1) ODE E-correction power n - (d-1) q_O = -sigma", sp.simplify(n - (d - 1) * qO + sig) == 0)
    check("(C1) Gamma E-correction power n - (D-1) q_G = sigma (p - 1) < 0 for d > 1",
          sp.simplify(n - (D - 1) * qG - sig * (p - 1)) == 0)
    check("(C1) q_O = p sigma/(1 - p)", sp.simplify(qO - p * sig / (1 - p)) == 0)
    ts = x ** qO / sig
    sad = sp.powdenest(sp.expand_power_base(p * x ** qG * sig ** p * ts ** (p - 1) / qG, force=True), force=True)
    check("(C2) tau* = x^(q_O)/sigma solves p x^(q_G) sigma^p tau^(p-1)/q_G = 1", sp.simplify(sad - 1) == 0, f"{sp.simplify(sad)}")
    ph = ts * v - x ** qG * (sig * ts * v) ** p / qG
    diff = sp.powdenest(sp.expand_power_base(sp.expand(ph - ts * (v - v ** p / p)), force=True), force=True)
    # run 3: sympy keeps x^(e1) x^(e2) ... as separate factors with uncancelled exponents (runs 1, 2 FAIL); rebuild each term with
    # the summed, cancel()-ed exponent of x and of v
    def norm(term):
        ex = {x: 0, v: 0}
        rest = sp.Integer(1)
        for f in sp.Mul.make_args(term):
            if f.is_Pow and f.base in ex:
                ex[f.base] += f.exp
            elif f in ex:
                ex[f] += 1
            else:
                rest *= f
        return rest * x ** sp.cancel(ex[x]) * v ** sp.cancel(ex[v])
    diff = sp.Add(*[norm(tm) for tm in sp.Add.make_args(sp.expand(diff))])
    check("(C2) tau - X^(q_G)/q_G = tau* (v - v^p/p) at tau = tau* v", sp.simplify(diff) == 0, f"{sp.simplify(diff)}")
    check("(C2) tau* h(1) = -x^(q_O)/q_O", sp.simplify(ts * (1 - 1 / p) + x ** qO / qO) == 0)
    hv = v - v ** p / p
    check("(C2) h'(1) = 0, h''(1) = 1 - p", sp.simplify(sp.diff(hv, v).subs(v, 1)) == 0
          and sp.simplify(sp.diff(hv, v, 2).subs(v, 1) - (1 - p)) == 0)


def c3():
    y = sp.Symbol("y", real=True)
    P = sp.Symbol("p", positive=True)
    t = sp.atan(y)
    R = sp.cos(P * t) / sp.cos(t) ** P - 1
    ser = sp.series(R, y, 0, 5).removeO()
    check("(C3) R(y) = p(1-p) y^2/2 + O(y^4)", sp.simplify(ser.coeff(y, 2) - P * (1 - P) / 2) == 0 and sp.simplify(ser.coeff(y, 1)) == 0
          and sp.simplify(ser.coeff(y, 3)) == 0, f"{sp.factor(ser.coeff(y, 2))}")
    # identity Re h(1+iy) - h(1) = -R/p, numerically at dps 40 (h(1) = 1 - 1/p)
    mp.mp.dps = 40
    worst = mp.mpf(0)
    for pp in (mp.mpf(2) / 5, mp.mpf(4) / 5):
        for yy in (mp.mpf("0.01"), mp.mpf("0.7"), mp.mpf(3), mp.mpf(1000)):
            vv = 1 + 1j * yy
            lhs = mp.re(vv - vv ** pp / pp) - (1 - 1 / pp)
            th = mp.atan(yy)
            rhs = -(mp.cos(pp * th) / mp.cos(th) ** pp - 1) / pp
            worst = max(worst, abs(lhs - rhs))
    check("(C3) Re h(1+iy) - h(1) = -R(y)/p (dps 40, two p, four y)", worst < mp.mpf(10) ** -35, mp.nstr(worst, 3))
    mp.mp.dps = 30
    for pr in (sp.Rational(2, 5), sp.Rational(1, 2), sp.Rational(4, 7), sp.Rational(2, 3), sp.Rational(3, 4), sp.Rational(4, 5)):
        pp = mp.mpf(pr.p) / pr.q
        mn = mp.inf
        for k in range(-60, 121):
            yy = mp.power(10, mp.mpf(k) / 20)
            th = mp.atan(yy)
            Rv = mp.cos(pp * th) / mp.cos(th) ** pp - 1
            mn = min(mn, Rv / min(yy ** 2, yy ** pp))
        print(f"   info (C3) p = {pr}: min over y in [1e-3, 1e6] of R/min(y^2, y^p) = {mp.nstr(mn, 6)}"
              f"  (p(1-p)/2 = {mp.nstr(pp * (1 - pp) / 2, 6)}, cos(p pi/2) = {mp.nstr(mp.cos(pp * mp.pi / 2), 6)})", flush=True)
        check(f"(C3) R/min(y^2, y^p) bounded below by a positive constant on the grid, p = {pr}", mn > mp.mpf("0.01"))


def c4():
    mp.mp.dps = 30
    pp = mp.mpf(4) / 5
    kap = mp.mpc("0.3", "-0.2")
    target = 1 / mp.sqrt(2 * mp.pi * (1 - pp))
    errs = []
    for lam in (50, 200, 800):
        f = lambda yy: mp.exp(lam * ((1 + 1j * yy) - (1 + 1j * yy) ** pp / pp - (1 - 1 / pp))) * (1 + 1j * yy) ** kap  # noqa: E731
        Iv = mp.quad(f, [-mp.inf, -1, -mp.mpf("0.1"), 0, mp.mpf("0.1"), 1, mp.inf])
        val = mp.sqrt(lam) / (2 * mp.pi) * Iv
        errs.append(abs(val / target - 1))
        print(f"   info (C4) lam = {lam}: sqrt(lam)/(2 pi) int = {mp.nstr(val, 12)}; |ratio - 1| = {mp.nstr(errs[-1], 4)}", flush=True)
    check("(C4) Laplace main term -> (2 pi (1-p))^(-1/2) with error ~ 1/lam (ratios 4 +- 1 between lam steps)",
          all(mp.mpf(3) < errs[i] / errs[i + 1] < mp.mpf(5) for i in range(2)), f"{[mp.nstr(e, 3) for e in errs]}")


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] != "S4":
        raise SystemExit("usage: spec1_s4_checks.py S4")
    SL.banner(__file__)
    c1_c2()
    c3()
    c4()
    SL.imported("imported(end)")
    print("S4CHK DONE: " + ("ALL PASS" if not FAIL else f"FAILURES {FAIL}"), flush=True)
