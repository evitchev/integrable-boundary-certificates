"""KIT HELPER (kit v1.2.6; not an archived file).  Quick recomputation of the SPEC1 S3 numerics at a few energies.

For each of the five points of the paper (Sol 3 a = 3, Sol 1 a = 3, Sol 3 a = 1, Sol 1 a = 1, Sol 1 a = 2) this calls the archived
solver path of spec1_s3.py (qpair: the byte-copied NUM2 ODE solver deps_snapshot/ode.py on both operators, at dps 70) at
E = 7 for momentum A and at E = 500 for momentum B, and checks:
  - the recomputed Q^O_j(E), Q^Gamma_j(-E) equal the archived out/s3/*.json values to >= 60 digits;
  - |omega Q^O_j/(m(theta_j) Q^Gamma_j) - 1| < 1e-60 with the closed-form omega of Theorem thm:spectral;
  - the tampers omega * 8/7 and m(theta_j) interchanged give deviations > 1e-3.
The full eight-energy runs with the dps-140 repeat are the full-mode checks spec1_s3_*.
Usage: spot_s3.py"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
STAGE = os.path.join(HERE, os.pardir, "spec1")
sys.path.insert(0, os.path.join(STAGE, "code"))
import spec1_s3 as S3  # noqa: E402
from spec1_s3 import SL, sp, mp  # noqa: E402

POINTS = [("Sol3", 3), ("Sol1", 3), ("Sol3", 1), ("Sol1", 1), ("Sol1", 2)]
CASES = [("A", 7), ("B", 500)]
FAIL = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""), flush=True)
    if not ok:
        FAIL.append(name)


def main():
    for sol, ai in POINTS:
        arch = json.load(open(os.path.join(STAGE, "out", "s3", f"{sol}_a{ai}.json")))
        for nm, E in CASES:
            mp.mp.dps = 70
            l0, l1 = S3.MOM[nm]
            pt, th, r, QO, QG, QF, _ = S3.qpair(sol, ai, l0, l1, E)
            rec = arch[nm][str(E)]
            aO = [mp.mpf(q) for q in rec["QO"]]
            aG = [mp.mpf(q) for q in rec["QG"]]
            nj = len(th)
            e_arch = max(max(abs(QO[j] / aO[j] - 1) for j in range(nj)), max(abs(QG[j] / aG[j] - 1) for j in range(nj)))
            w = mp.mpf(sp.N(SL.omega_sealed(pt, r), 90))
            ms = [mp.mpf(sp.N(SL.m_sym(pt, r, x), 90)) for x in th]
            perm = list(range(1, nj)) + [0]
            dev = max(abs(w * QO[j] / (ms[j] * QG[j]) - 1) for j in range(nj))
            tw = min(abs(w * mp.mpf(8) / 7 * QO[j] / (ms[j] * QG[j]) - 1) for j in range(nj))
            tp = max(abs(w * QO[j] / (ms[perm[j]] * QG[j]) - 1) for j in range(nj))
            tag = f"{sol} a={ai} t={pt['t']} {nm} E={E}"
            check(f"{tag}: recomputed Q's = archived out/s3 values (>= 60 digits)", e_arch < mp.mpf(10) ** -60)
            check(f"{tag}: relation |omega Q^O/(m Q^Gamma) - 1| < 1e-60", dev < mp.mpf(10) ** -60)
            check(f"{tag}: tampers omega*8/7 and m interchanged fail (> 1e-3)", tw > mp.mpf(10) ** -3 and tp > mp.mpf(10) ** -3,
                  f"{mp.nstr(tw, 2)}, {mp.nstr(tp, 2)}")
    print("SPOT S3: " + ("ALL PASS" if not FAIL else f"FAILURES {FAIL}"), flush=True)
    return 0 if not FAIL else 1


if __name__ == "__main__":
    sys.exit(main())
