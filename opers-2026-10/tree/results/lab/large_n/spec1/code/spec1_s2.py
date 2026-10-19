"""SPEC1 S2 (SPEC1_PREREG sec. 4 sealing procedure): checks of the SEALED omega = sigma^(g-1/2) sqrt((d+a)/(2 pi (d-1))).  NEW FILE.
Usage: spec1_s2.py S2
S2a (symbolic/exact):
  - the power identity  -g q_O + beta_G (1 + q_O/sigma) + q_O/2 = beta_O, symbolic in (a, d, s, r, sigma);
  - the saddle: f'(tau*) = 0, f(tau*) = -x^q_O/q_O, f''(tau*) = (1-p)/tau* at tau* = x^q_O/sigma (x > 0), symbolic;
  - beta_G, beta_O by EXACT formal conjugation of each operator (x^beta e^(-u/q), u = x^q; theta u^j = q j u^j), at every S0 point
    a = 1..6 of both solutions, symbolic in l0, l1, compared with the PREREG formulas.  Tamper: beta + 1/7 must leave the
    u^(D-1) (resp. u^(d-1)) coefficient non-zero.
S2b: omega(t = 2) = 3^(-1/6) pi^(-1/2) exactly.
S2c (independent route, E = 0): Meijer-G connection coefficients of both forms (spec1_lib), validated against the archived
  t = 2 closed forms (PROOF3/NUM2 ops.exact_three_term_E0, exact_sixth_E0; byte copy) at NUM2's points A, B; then at every S0 point
  a = 1..6 (both solutions) and momenta A, B: omega * Q^O_h(0) / (m(theta_h) Q^G_h(0)) = 1 for every dynamical h, at 60 digits
  (and 120 digits at the S3 points).  Exponent consistency: c*vt = beta (exact rationals).  Tamper: omega from g + 1/7 must fail."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import spec1_lib as SL  # noqa: E402
from spec1_lib import sp, mp, F  # noqa: E402
import ops as NUM2OPS  # noqa: E402  (byte copy of PROOF3/NUM2 ops.py)

FAIL = []
MOM = {"A": (sp.Rational(1, 3), sp.Rational(3, 5)), "B": (sp.Rational(9, 10), sp.Rational(1, 5))}
S3PTS = [("Sol3", 3), ("Sol1", 3), ("Sol3", 1), ("Sol1", 1), ("Sol1", 2), ("Sol3", 2)]


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""), flush=True)
    if not ok:
        FAIL.append(name)


def s2a_symbolic():
    a, d, s, r, sig, x = sp.symbols("a d s r sigma x", positive=True)
    c = sig * (a + 1)
    D = d + a
    qG, qO = c / D, c / (d - 1)
    p = qG / sig
    g = (s - r) / sig + (1 - a) / 2
    bG = s - qG * (D - 1) / 2 - a * r / D
    bO = s - c / 2 + r / (d - 1)
    ident = sp.simplify(-g * qO + bG * (1 + qO / sig) + qO / 2 - bO)
    check("S2a power identity -g q_O + beta_G (1 + q_O/sigma) + q_O/2 = beta_O (symbolic)", ident == 0, f"residual {ident}")
    tau = sp.Symbol("tau", positive=True)
    f = tau - x ** qG * (sig * tau) ** p / qG
    ts = x ** qO / sig

    def red(e):
        """run 2: denest (x^A)^B (x > 0) before simplifying (run 1: sympy left the nested powers, POWDENEST)."""
        return sp.simplify(sp.powsimp(sp.powdenest(sp.expand_power_base(e, force=True), force=True), force=True))
    fp = red(sp.diff(f, tau).subs(tau, ts))
    fv = red(f.subs(tau, ts) + x ** qO / qO)
    fpp = red(sp.diff(f, tau, 2).subs(tau, ts) - (1 - p) / ts)
    check("S2a saddle: f'(tau*) = 0, f(tau*) = -x^q_O/q_O, f''(tau*) = (1-p)/tau* (symbolic, x > 0)",
          fp == 0 and fv == 0 and fpp == 0, f"{fp}, {fv}, {fpp}")
    Xs = red(x * (sig * ts) ** (1 / sig))
    check("S2a X* = x^(1 + q_O/sigma)", red(Xs / x ** (1 + qO / sig)) == 1, f"{Xs}")


def conj_coeff(roots, q, sb, D, extra):
    """Pi(T') 1 with T' = sb + u - theta (theta u^j = q j u^j), roots of Pi in T; returns {j: coeff}; then 'extra' (a dict) is
    subtracted."""
    Fd = {0: sp.Integer(1)}
    for rho in roots:
        out = {}
        for j, f in Fd.items():
            out[j] = out.get(j, 0) + f * (sb - rho - q * j)
            out[j + 1] = out.get(j + 1, 0) + f
        Fd = out
    for j, v in extra.items():
        Fd[j] = Fd.get(j, 0) - v
    return {j: sp.expand(v) for j, v in Fd.items()}


def s2a_beta():
    l0, l1, beta = sp.symbols("l0 l1 beta")
    ok = okt = True
    for sol in ("Sol3", "Sol1"):
        for ai in range(1, 7):
            pt = SL.point(sol, ai)
            s, sig, c = (sp.Rational(pt[k]) for k in ("s", "sigma", "c"))
            d, D = pt["d"], pt["D"]
            if sol == "Sol3":
                Proots, r = [l0, -l0, l1, -l1], sp.Integer(0)
            else:
                Proots, r = [l1, -l1], l0
            Groots = [r + sig * sp.Rational(ai - 1, 2) - sig * j for j in range(ai)]
            # Gamma form: Pi(T')1 - x^c, x^c = u^D (q = c/D)
            qG = c / D
            cg = conj_coeff(Proots + Groots, qG, s - beta, D, {D: 1})
            assert sp.expand(cg.get(D, 0)) == 0
            bG = sp.solve(cg[D - 1], beta)
            # ODE: P(T')1 - u^(d-1) (T' - r - c/2) 1 = P(T')1 - u^d - (s - beta - r - c/2) u^(d-1)
            qO = c / (d - 1)
            co = conj_coeff(Proots, qO, s - beta, d, {d: 1, d - 1: s - beta - r - c / 2})
            assert sp.expand(co.get(d, 0)) == 0
            bO = sp.solve(co[d - 1], beta)
            fG, fO = SL.beta_gamma(pt, r), SL.beta_ode(pt, r)
            okp = len(bG) == 1 and len(bO) == 1 and sp.simplify(bG[0] - fG) == 0 and sp.simplify(bO[0] - fO) == 0
            ok &= okp
            okt &= sp.expand(cg[D - 1].subs(beta, fG + sp.Rational(1, 7))) != 0 and \
                sp.expand(co[d - 1].subs(beta, fO + sp.Rational(1, 7))) != 0
            print(f"   {sol} a={ai}: beta_G = {bG} (formula {fG}); beta_O = {bO} (formula {fO})", flush=True)
    check("S2a beta_G, beta_O by exact conjugation = PREREG formulas, a = 1..6, both solutions (symbolic l0, l1)", ok)
    check("S2a tamper beta + 1/7 leaves the matching coefficient non-zero", okt)


def s2b():
    pt = SL.point("Sol3", 2)
    w = SL.omega_sealed(pt, sp.Integer(0))
    target = sp.Integer(3) ** sp.Rational(-1, 6) / sp.sqrt(sp.pi)
    check("S2b omega(t = 2) = 3^(-1/6) pi^(-1/2) exactly", sp.simplify(w / target) == 1, f"{w}")


def s2c():
    mp.mp.dps = 60
    print("== S2c positive control: generic Meijer formulas vs archived t = 2 closed forms (ops.py byte copy)", flush=True)
    pt = SL.point("Sol3", 2)
    okc = True
    for nm, (l0, l1) in MOM.items():
        th, r, frozen = SL.exponents(pt, l0, l1)
        qo, _ = SL.meijer_Q0_ode(pt, th, r)
        qg, _ = SL.meijer_Q0_gamma(pt, th, frozen)
        ao, _ = NUM2OPS.exact_three_term_E0(2, mp.mpf(l0.p) / l0.q, mp.mpf(l1.p) / l1.q)
        ag, _ = NUM2OPS.exact_sixth_E0(mp.mpf(l0.p) / l0.q, mp.mpf(l1.p) / l1.q)
        # archived orders: three_term [s-l0, s+l0, s-l1, s+l1]; sixth [1, 4, s-l0, s+l0, s-l1, s+l1]
        eo = max(abs(qo[h] / ao[h] - 1) for h in range(4))
        eg = max(abs(qg[h] / ag[h + 2] - 1) for h in range(4))
        okc &= eo < mp.mpf(10) ** -55 and eg < mp.mpf(10) ** -55
        print(f"   t=2 point {nm}: max rel diff ODE {mp.nstr(eo, 3)}, Gamma {mp.nstr(eg, 3)}", flush=True)
    check("S2c control: generic E = 0 formulas reproduce the archived t = 2 closed forms (both forms, A and B)", okc)
    print("== S2c ratio test of the SEALED omega", flush=True)
    okall = oktam = okexp = True
    worst = mp.mpf(0)
    for sol in ("Sol3", "Sol1"):
        for ai in range(1, 7):
            pt = SL.point(sol, ai)
            for nm, (l0, l1) in MOM.items():
                th, r, frozen = SL.exponents(pt, l0, l1)
                # exact exponent consistency c*vt = beta
                c = sp.Rational(pt["c"])
                d, D = pt["d"], pt["D"]
                vtO = (sp.Rational(1 - (d - 1), 2) + sum(x / c for x in th) - (sp.Rational(pt["s"]) + c / 2 - r) / c) / (d - 1)
                vtG = (sp.Rational(1 - D, 2) + sum(x / c for x in list(th) + list(frozen))) / D
                okexp &= sp.simplify(c * vtO - SL.beta_ode(pt, r)) == 0 and sp.simplify(c * vtG - SL.beta_gamma(pt, r)) == 0
                for dps in ((60, 120) if (sol, ai) in S3PTS else (60,)):
                    mp.mp.dps = dps
                    qo, _ = SL.meijer_Q0_ode(pt, th, r)
                    qg, _ = SL.meijer_Q0_gamma(pt, th, frozen)
                    w = mp.mpf(sp.N(SL.omega_sealed(pt, r), dps + 10))
                    sig = sp.Rational(pt["sigma"])
                    gt = SL.g_of(pt, r) + sp.Rational(1, 7)
                    wt = mp.mpf(sp.N(sig ** (gt - sp.Rational(1, 2)) * sp.sqrt(sp.Rational(d + ai, 2 * (d - 1)) / sp.pi), dps + 10))
                    devs, devt = [], []
                    for h, thh in enumerate(th):
                        mth = mp.mpf(sp.N(SL.m_sym(pt, r, thh), dps + 10))
                        R = w * qo[h] / (mth * qg[h])
                        devs.append(abs(R - 1))
                        devt.append(abs(wt * qo[h] / (mth * qg[h]) - 1))
                    dev = max(devs)
                    tol = mp.mpf(10) ** -(dps - 10)
                    okall &= dev < tol
                    oktam &= min(devt) > mp.mpf(10) ** -3
                    worst = max(worst, dev * mp.mpf(10) ** (dps - 60))
                    print(f"   {sol} a={ai} t={pt['t']} {nm} dps {dps}: max_h |omega Q^O/(m Q^G) - 1| = {mp.nstr(dev, 3)}; "
                          f"tamper min dev {mp.nstr(min(devt), 3)}; 1/omega = {mp.nstr(1 / w, 25)}", flush=True)
    check("S2c exact exponent consistency c*vt = beta (both forms, every point/momentum)", okexp)
    check("S2c sealed omega: omega Q^O_h(0) = m(theta_h) Q^G_h(0) for every h, every S0 point a = 1..6, both solutions, A and B "
          "(< 10^-(dps-10))", okall)
    check("S2c tamper (omega from g + 1/7) fails everywhere (> 1e-3)", oktam)


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] != "S2":
        raise SystemExit("usage: spec1_s2.py S2")
    SL.banner(__file__)
    s2a_symbolic()
    s2a_beta()
    s2b()
    s2c()
    SL.imported("imported(end)")
    print("S2 DONE: " + ("ALL PASS" if not FAIL else f"FAILURES {FAIL}"), flush=True)
