"""SPEC1 S3 (SPEC1_PREREG sec. 5 + lead amendment): numerical test of Q^O_j(E) = (m(theta_j)/omega) Q^G_j(E_Gamma = -E) with the
SEALED omega.  NEW FILE.
Usage: spec1_s3.py S3 <Sol3|Sol1> <a> [ctrl]
  ctrl (only Sol3 2): the positive control against PROOF3/NUM2's archived t1.json (Q3, Q6 at E in {-20, -5, 2, 50, 500}).
Operators (byte-copied NUM2 Op solver; equations made monic in theta by the factor (-1)^order):
  L_O = P(T) - x^c (T - r - c/2) - E x^n;  L_G = P(T) G_a(T - r) - x^c - E x^n  (same E; i.e. the Gamma form at E_Gamma = -E).
Energies E in {-20, -5, 0, 1, 7, 37, 100, 500}, momenta A (1/3, 3/5), B (9/10, 1/5); mp.dps = 70; E = 37 repeated at dps 140.
PASS: |omega Q^O_j/(m(theta_j) Q^G_j) - 1| < 1e-30 for every j, E, momentum.  E = 0: both Q's against the Meijer closed forms
(>= 40 digits).  Tampers: T-omega (omega (1 + 1/7)), T-perm (m(theta) with exponents permuted: swap for Sol 1, cyclic shift for
Sol 3) must fail; T-order (ordering x^c T in L_O) is information."""
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import spec1_lib as SL  # noqa: E402
from spec1_lib import sp, mp, F  # noqa: E402
from ode import Op  # noqa: E402  (byte copy of NUM2/PROOF3 ode.py)

FAIL = []
MOM = {"A": (sp.Rational(1, 3), sp.Rational(3, 5)), "B": (sp.Rational(9, 10), sp.Rational(1, 5))}
ES = [-20, -5, 0, 1, 7, 37, 100, 500]
SPECIAL = {("Sol3", 1): "SPECIAL (t = 5/3: Sol 2 = Sol 3 locus)", ("Sol1", 1): "SPECIAL (t = 3, N = -2 anchor)",
           ("Sol1", 2): "SPECIAL (t = 5, N = 0)", ("Sol3", 2): "CONTROL (t = 2, PROOF3 point)"}
OUT = os.path.join(SL.ROOT, "out", "s3")


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""), flush=True)
    if not ok:
        FAIL.append(name)


def mpq(x):
    x = sp.Rational(x)
    return mp.mpf(x.p) / x.q


def theta_poly(roots_T, s):
    """prod_i (T - rho_i) with T = s - theta, times (-1)^deg so that it is monic in theta: prod_i (theta - (s - rho_i)).
    Returns ascending coefficient list (mpf) and the list of theta-roots."""
    th_roots = [s - rho for rho in roots_T]
    P = [mp.mpf(1)]
    for e in th_roots:
        ev = mp.mpf(sp.N(e, mp.mp.dps + 20))
        new = [mp.mpf(0)] * (len(P) + 1)
        for i, cf in enumerate(P):
            new[i + 1] += cf
            new[i] -= ev * cf
        P = new
    return P, th_roots


def build(sol, ai, l0, l1, E, wrong_order=False):
    pt = SL.point(sol, ai)
    s, sig, c, n = (sp.Rational(pt[k]) for k in ("s", "sigma", "c", "n"))
    d, D = pt["d"], pt["D"]
    th, r, frozen = SL.exponents(pt, l0, l1)
    Proots = [s - x for x in th]                                  # T-roots of P
    Groots = [r + sig * sp.Rational(ai - 1, 2) - sig * j for j in range(ai)]
    E = mp.mpf(E)
    # ODE, monic in theta: (-1)^d [P(T) - x^c (s - theta - r - c/2) - E x^n]
    Po, _ = theta_poly(Proots, s)
    sd = (-1) ** d
    shift = sp.Integer(0) if wrong_order else c / 2
    R1 = [-sd * mpq(0) + sd * (-(mp.mpf(sp.N(s - r - shift, mp.mp.dps + 20)))), sd * mp.mpf(1)]   # (-1)^d (theta - (s - r - shift))
    opO = Op(Po, [(F(pt["c"]), R1), (F(pt["n"]), [-sd * E])])
    # Gamma form, monic in theta: (-1)^D [Pi(T) - x^c - E x^n]
    Pg, gth = theta_poly(Proots + Groots, s)
    sD = (-1) ** D
    opG = Op(Pg, [(F(pt["c"]), [-sD * mp.mpf(1)]), (F(pt["n"]), [-sD * E])])
    thO = [mp.mpf(sp.N(x, mp.mp.dps + 20)) for x in th]
    thG = [mp.mpf(sp.N(x, mp.mp.dps + 20)) for x in list(th) + list(frozen)]
    return pt, opO, thO, opG, thG, th, r


def qpair(sol, ai, l0, l1, E, wrong_order=False):
    pt, opO, thO, opG, thG, th, r = build(sol, ai, l0, l1, E, wrong_order)
    x0 = 2 / (abs(mp.mpf(E)) + 1) ** (1 / mpq(pt["n"]))
    t0 = time.time()
    _, QO = opO.Q(thO, x0)
    t1 = time.time()
    _, QG = opG.Q(thG, x0)
    t2 = time.time()
    return pt, th, r, QO, QG[:len(thO)], QG[len(thO):], (t1 - t0, t2 - t1)


def run_point(sol, ai):
    os.makedirs(OUT, exist_ok=True)
    lab = SPECIAL.get((sol, ai), "GENERIC")
    pt = SL.point(sol, ai)
    print(f"== {sol} a={ai} t={pt['t']} [{lab}]: ODE order {pt['d']}, Gamma order {pt['D']}, sigma {pt['sigma']}, n {pt['n']}, "
          f"c {pt['c']}", flush=True)
    res = {}
    worst = mp.mpf(0)
    tam_w = tam_p = mp.inf
    order_info = []
    for nm, (l0, l1) in MOM.items():
        res[nm] = {}
        for E in ES:
            mp.mp.dps = 70
            pt, th, r, QO, QG, QF, (tO, tG) = qpair(sol, ai, l0, l1, E)
            w = mp.mpf(sp.N(SL.omega_sealed(pt, r), 90))
            ms = [mp.mpf(sp.N(SL.m_sym(pt, r, x), 90)) for x in th]
            dev = [abs(w * QO[j] / (ms[j] * QG[j]) - 1) for j in range(len(th))]
            perm = list(range(1, len(th))) + [0]
            dT_w = [abs(w * (1 + mp.mpf(1) / 7) * QO[j] / (ms[j] * QG[j]) - 1) for j in range(len(th))]
            dT_p = [abs(w * QO[j] / (ms[perm[j]] * QG[j]) - 1) for j in range(len(th))]
            worst = max(worst, max(dev))
            tam_w, tam_p = min(tam_w, min(dT_w)), min(tam_p, max(dT_p))
            res[nm][str(E)] = {"QO": [str(q) for q in QO], "QG": [str(q) for q in QG], "QG_frozen": [str(q) for q in QF],
                               "dev": [str(x) for x in dev], "seconds": [tO, tG]}
            print(f"   {nm} E={E}: max_j |omega Q^O/(m Q^G) - 1| = {mp.nstr(max(dev), 3)}; T-omega min {mp.nstr(min(dT_w), 3)}; "
                  f"T-perm max {mp.nstr(max(dT_p), 3)}; time ODE {tO:.1f}s Gamma {tG:.1f}s", flush=True)
            if E == 0:
                qo0, _ = SL.meijer_Q0_ode(pt, th, r)
                _, rr, frozen = SL.exponents(pt, l0, l1)
                qg0, _ = SL.meijer_Q0_gamma(pt, th, frozen)
                e0 = max(max(abs(QO[j] / qo0[j] - 1) for j in range(len(th))), max(abs(QG[j] / qg0[j] - 1) for j in range(len(th))))
                check(f"S3 E=0 control {sol} a={ai} {nm}: both solvers = Meijer closed forms (>= 40 digits)", e0 < mp.mpf(10) ** -40,
                      f"max rel diff {mp.nstr(e0, 3)}")
            if E == 37:
                mp.mp.dps = 140
                pt2, th2, r2, QO2, QG2, _, (tO2, tG2) = qpair(sol, ai, l0, l1, E)
                e2 = max(max(abs(QO2[j] / QO[j] - 1) for j in range(len(th))), max(abs(QG2[j] / QG[j] - 1) for j in range(len(th))))
                w2 = mp.mpf(sp.N(SL.omega_sealed(pt2, r2), 160))
                ms2 = [mp.mpf(sp.N(SL.m_sym(pt2, r2, x), 160)) for x in th2]
                dev2 = max(abs(w2 * QO2[j] / (ms2[j] * QG2[j]) - 1) for j in range(len(th2)))
                check(f"S3 dps-140 repeat {sol} a={ai} {nm} E=37: Q's agree with dps 70 to >= 50 digits", e2 < mp.mpf(10) ** -50,
                      f"{mp.nstr(e2, 3)}; dps-140 relation deviation {mp.nstr(dev2, 3)} ({tO2:.0f}s + {tG2:.0f}s)")
                mp.mp.dps = 70
                # T-order (information)
                _, thw, rw, QOw, QGw, _, _ = qpair(sol, ai, l0, l1, E, wrong_order=True)
                order_info.append(max(abs(w * QOw[j] / (ms[j] * QGw[j]) - 1) for j in range(len(thw))))
    check(f"S3 {sol} a={ai} [{lab}]: relation holds, max deviation < 1e-30 over E, j, A/B", worst < mp.mpf(10) ** -30,
          f"worst {mp.nstr(worst, 3)}")
    check(f"S3 {sol} a={ai}: T-omega fails (min deviation > 1e-3)", tam_w > mp.mpf(10) ** -3, f"{mp.nstr(tam_w, 3)}")
    check(f"S3 {sol} a={ai}: T-perm fails (deviation > 1e-3 at every E)", tam_p > mp.mpf(10) ** -3, f"min over E {mp.nstr(tam_p, 3)}")
    print(f"   info T-order (x^c T ordering) at E = 37: deviations {[mp.nstr(x, 3) for x in order_info]}", flush=True)
    json.dump(res, open(os.path.join(OUT, f"{sol}_a{ai}.json"), "w"), indent=1)


def run_control():
    mp.mp.dps = 60
    t1 = json.load(open(os.path.join(SL.ROOT, "inputs", "proof3", "t1.json")))
    worst = mp.mpf(0)
    for nm in ("A", "B"):
        l0, l1 = MOM[nm]
        for E in (-20, -5, 2, 50, 500):
            pt, th, r, QO, QG, QF, (tO, tG) = qpair("Sol3", 2, l0, l1, E)
            arch = t1[nm][str(mp.mpf(E))]
            aO = [mp.mpf(x) for x in arch["Q3"]]
            aG = [mp.mpf(x) for x in arch["Q6"]]       # archived order: frozen (1, 4) then dynamical
            e = max(max(abs(QO[j] / aO[j] - 1) for j in range(4)), max(abs(QG[j] / aG[j + 2] - 1) for j in range(4)),
                    max(abs(QF[k] / aG[k] - 1) for k in range(2)))
            worst = max(worst, e)
            print(f"   control t=2 {nm} E={E}: max rel diff vs t1.json (Q3, Q6 incl. frozen) {mp.nstr(e, 3)} ({tO:.1f}s + {tG:.1f}s)",
                  flush=True)
    check("S3 positive control: Sol 3 a = 2 reproduces the archived t1.json Q3/Q6 (>= 45 digits)", worst < mp.mpf(10) ** -45,
          f"worst {mp.nstr(worst, 3)}")


if __name__ == "__main__":
    if len(sys.argv) < 4 or sys.argv[1] != "S3":
        raise SystemExit("usage: spec1_s3.py S3 <Sol3|Sol1> <a> [ctrl]")
    SL.banner(__file__)
    sol, ai = sys.argv[2], int(sys.argv[3])
    T0 = time.time()
    if len(sys.argv) > 4 and sys.argv[4] == "ctrl":
        assert (sol, ai) == ("Sol3", 2)
        run_control()
    else:
        run_point(sol, ai)
    print(f"total wall time {time.time() - T0:.0f}s", flush=True)
    SL.imported("imported(end)")
    print("S3 DONE: " + ("ALL PASS" if not FAIL else f"FAILURES {FAIL}"), flush=True)
