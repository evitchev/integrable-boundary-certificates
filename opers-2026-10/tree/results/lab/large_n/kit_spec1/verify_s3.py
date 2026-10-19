"""KIT HELPER (kit v1.2.6; not an archived file).  Fast verifier of the archived SPEC1 S3 outputs (results/lab/large_n/spec1/out/s3/*.json).

The archived S3 runs (spec1_s3.py, 0.5-2 h per point, dominated by the dps-140 repeat at E = 37) stored, for each of the five points,
both momenta (A = (1/3, 3/5), B = (9/10, 1/5)) and E in {-20, -5, 0, 1, 7, 37, 100, 500}, the dps-70 Q-functions Q^O_j(E) of the
three-term ODE and Q^Gamma_j(-E) of the Gamma form.  This script, without re-solving any ODE:
  - checks the sha256 of each JSON against ORIGIN_SHA256 (the stage's sealed SHA256SUMS);
  - recomputes omega (closed form of Theorem thm:spectral) and m(theta_j) exactly from spec1_lib.py and the deviation
    |omega Q^O_j/(m(theta_j) Q^Gamma_j) - 1| at dps 70, and compares it with the archived 'dev' field;
  - checks the paper's statements: deviation < 1e-30 everywhere (the registered criterion) and < 1e-56 everywhere, < 1e-64 except
    at E = 0 for Solution 1; at E = 0 the Q's against the exact Meijer-G values (spec1_lib), whose difference matches the
    relation's deviation;
  - the two tampers of the paper: omega * 8/7, and the multipliers m(theta_j) interchanged between exponents (cyclic shift for
    Sol 3, swap for Sol 1, as in spec1_s3.py); each must give a deviation > 1e-3 at every energy.
Not covered here: the dps-140 repeat at E = 37 (only its summary is in the logs); the full re-runs (s3_* checks, full mode) cover it.
Usage: verify_s3.py"""
import hashlib
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
STAGE = os.path.join(HERE, os.pardir, "spec1")
sys.path.insert(0, os.path.join(STAGE, "code"))
import spec1_lib as SL  # noqa: E402
from spec1_lib import sp, mp  # noqa: E402

MOM = {"A": (sp.Rational(1, 3), sp.Rational(3, 5)), "B": (sp.Rational(9, 10), sp.Rational(1, 5))}
ES = [-20, -5, 0, 1, 7, 37, 100, 500]
POINTS = [("Sol3", 3), ("Sol1", 3), ("Sol3", 1), ("Sol1", 1), ("Sol1", 2)]
FAIL = []
E0SOL1 = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""), flush=True)
    if not ok:
        FAIL.append(name)


def e3(x):
    return mp.nstr(x, 2)


def main():
    want = {}
    for line in open(os.path.join(HERE, "ORIGIN_SHA256")):
        if line.strip():
            h, f = line.split()
            want[os.path.basename(f)] = h
    for sol, ai in POINTS:
        fn = f"{sol}_a{ai}.json"
        path = os.path.join(STAGE, "out", "s3", fn)
        h = hashlib.sha256(open(path, "rb").read()).hexdigest()
        check(f"{fn}: sha256 = ORIGIN_SHA256", h == want.get(fn), h[:16])
        d = json.load(open(path))
        pt = SL.point(sol, ai)
        print(f"== {sol} a={ai} t={pt['t']}: ODE order {pt['d']}, Gamma order {pt['D']}", flush=True)
        worst = worst_ne0 = mp.mpf(0)
        same_dev = 0
        n_dev = 0
        tw_min = tp_min = mp.inf
        for nm, (l0, l1) in MOM.items():
            mp.mp.dps = 70
            th, r, frozen = SL.exponents(pt, l0, l1)
            w = mp.mpf(sp.N(SL.omega_sealed(pt, r), 90))
            ms = [mp.mpf(sp.N(SL.m_sym(pt, r, x), 90)) for x in th]
            nj = len(th)
            perm = list(range(1, nj)) + [0]
            row = []
            for E in ES:
                rec = d[nm][str(E)]
                QO = [mp.mpf(q) for q in rec["QO"]]
                QG = [mp.mpf(q) for q in rec["QG"]]
                dev = [abs(w * QO[j] / (ms[j] * QG[j]) - 1) for j in range(nj)]
                for j in range(nj):
                    n_dev += 1
                    a_dev = mp.mpf(rec["dev"][j])
                    same_dev += (abs(dev[j] - a_dev) <= mp.mpf(10) ** -3 * a_dev + mp.mpf(10) ** -69)
                tw = min(abs(w * mp.mpf(8) / 7 * QO[j] / (ms[j] * QG[j]) - 1) for j in range(nj))
                tp = max(abs(w * QO[j] / (ms[perm[j]] * QG[j]) - 1) for j in range(nj))
                tw_min, tp_min = min(tw_min, tw), min(tp_min, tp)
                worst = max(worst, max(dev))
                if E != 0:
                    worst_ne0 = max(worst_ne0, max(dev))
                row.append(f"{E}:{e3(max(dev))}")
                if E == 0:
                    qo0, _ = SL.meijer_Q0_ode(pt, th, r)
                    qg0, _ = SL.meijer_Q0_gamma(pt, th, frozen)
                    eo = max(abs(QO[j] / qo0[j] - 1) for j in range(nj))
                    eg = max(abs(QG[j] / qg0[j] - 1) for j in range(nj))
                    rel0 = max(dev)
                    if sol == "Sol1":
                        E0SOL1.append(rel0)
                    print(f"   {nm} E=0: relation dev {e3(rel0)}; Q^O vs Meijer {e3(eo)}; Q^Gamma vs Meijer {e3(eg)}", flush=True)
                    check(f"{sol} a={ai} {nm} E=0: both Q's = exact Meijer-G values to >= 40 digits", max(eo, eg) < mp.mpf(10) ** -40)
                    check(f"{sol} a={ai} {nm} E=0: relation deviation is of the size of the Meijer difference (ratio in [1/10, 10])",
                          max(eo, eg) > 0 and mp.mpf(1) / 10 <= rel0 / max(eo, eg) <= 10 if rel0 > mp.mpf(10) ** -64 else True,
                          "(only tested where dev > 1e-64)")
            print(f"   {nm} max_j dev by E: " + " ".join(row), flush=True)
        print(f"   recomputed deviations = archived 'dev': {same_dev}/{n_dev}", flush=True)
        check(f"{sol} a={ai}: recomputed deviations equal the archived ones", same_dev == n_dev)
        check(f"{sol} a={ai}: relation holds, max deviation < 1e-30 (registered criterion)", worst < mp.mpf(10) ** -30, f"worst {e3(worst)}")
        check(f"{sol} a={ai}: max deviation < 2e-56 (paper: 'about 1e-56 or better')", worst < 2 * mp.mpf(10) ** -56, f"worst {e3(worst)}")
        if sol == "Sol1":
            check(f"{sol} a={ai}: max deviation at E != 0 < 1e-64 (paper: better than 1e-64 except at E = 0 for Solution 1)",
                  worst_ne0 < mp.mpf(10) ** -64, f"worst {e3(worst_ne0)}")
        else:
            check(f"{sol} a={ai}: max deviation at every E < 1e-64 (paper: better than 1e-64 except at E = 0 for Solution 1)",
                  worst < mp.mpf(10) ** -64, f"worst {e3(worst)}")
        check(f"{sol} a={ai}: tamper omega * 8/7 fails at every E (min dev > 1e-3)", tw_min > mp.mpf(10) ** -3, f"{e3(tw_min)} (= 1/7 by construction)")
        check(f"{sol} a={ai}: tamper m(theta_j) interchanged fails at every E (min over E of max_j dev > 1e-3)", tp_min > mp.mpf(10) ** -3,
              f"{e3(tp_min)}")
    lo, hi = min(E0SOL1), max(E0SOL1)
    check("Solution 1, E = 0 (3 points x 2 momenta): deviations between about 1e-61 and 1e-56 (paper; tested as [1e-63, 2e-56])",
          len(E0SOL1) == 6 and mp.mpf(10) ** -63 <= lo and hi <= 2 * mp.mpf(10) ** -56, f"range {e3(lo)} .. {e3(hi)}")
    print("VERIFY S3 ARCHIVE: " + ("ALL PASS" if not FAIL else f"FAILURES {FAIL}"), flush=True)
    return 0 if not FAIL else 1


if __name__ == "__main__":
    sys.exit(main())
