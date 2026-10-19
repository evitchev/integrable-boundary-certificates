"""SPEC1 S0 (scope), exact.  NEW FILE.  Usage: spec1_s0.py S0
Hypotheses (brief): positive integer a, sigma = n/a > 0, deg L = 1, deg P >= 2 (and no resonance at chosen momenta).
Parameters from the paper (opers_npb.tex, sec. "The opers"):
  p^2 = 2(t-1)/(t+1), k = 2(t-1)/(3-t), b as tabulated, a = b/(1-b), sigma = n/a, s = (n-1)/2, c = n/b.
  Sol 1: P = T^2 - l1^2, L = T - l0, n = (t+3)/2, b = p^2/2  (a = (t-1)/2).
  Sol 2: P = T(T^2 - l1^2), L = T^2 - l0^2 (deg L = 2: OUT of scope), n = 2k+3, b = p^2 (a = k).
  Sol 3: P = (T^2-l0^2)(T^2-l1^2), L = T, n = k+4, b = p^2 (a = k).
Derived per point (deg L = 1, root r): D = deg P + a (order of the integer-a Gamma form), g = (s - r)/sigma + (1-a)/2,
frozen lattice base_r = s - r + sigma (1-a)/2 (zeros of m: base + sigma Z>=0; frozen exponents base + sigma j, j < a),
p = (a+1)/(deg P + a) (Hankel convergence needs p < 1 <=> deg P >= 2), q_ODE = c/(deg P - 1), q_Gamma = c/D.
No-resonance (H) at a momentum point: theta_i - theta_j not in sigma Z (i != j), theta_j not in base + sigma Z (exact rationals)."""
import sys
from fractions import Fraction as F

DATA = {  # record numerical Q-function data (archives): (solution, t) -> source
    ("Sol3", F(2)): "NUM2 item 325 (results/lab/large_n/num2, k = 2; also the sixth-order Gamma form there, T1)",
    ("Sol3", F(9, 4)): "NUM2 item 325 (k = 10/3, non-integer a)",
    ("Sol2", F(2)): "NUM3 item 329 (k = 2; deg L = 2)",
    ("Sol2", F(9, 4)): "NUM3 item 329 (k = 10/3)",
    ("Sol1", F(2)): "NUM1 (M = 5, i.e. t = 2: a = 1/2, non-integer)",
    ("Sol1", F(81, 41)): "NUM1 (M = 51/10, i.e. t = 81/41: a = 20/41, non-integer)",  # run 2: was F(51, 41) (wrong t)
}


def params(sol, t):
    t = F(t)
    p2 = 2 * (t - 1) / (t + 1)
    if sol == "Sol1":
        n, b, degP, degL = (t + 3) / 2, p2 / 2, 2, 1
    elif sol == "Sol2":
        k = 2 * (t - 1) / (3 - t)
        n, b, degP, degL = 2 * k + 3, p2, 3, 2
    else:
        k = 2 * (t - 1) / (3 - t)
        n, b, degP, degL = k + 4, p2, 4, 1
    a = b / (1 - b)
    sigma = n / a
    s = (n - 1) / 2
    c = n / b
    return dict(t=t, n=n, b=b, a=a, sigma=sigma, s=s, c=c, degP=degP, degL=degL)


def t_of(sol, a):
    if sol == "Sol1":
        return 2 * F(a) + 1
    return (3 * F(a) + 2) / (F(a) + 2)          # k = a for Sol 2 and Sol 3


def check_H(sol, pr, l0, l1):
    s, sigma, a = pr["s"], pr["sigma"], pr["a"]
    if sol == "Sol3":
        r = F(0)
        th = [s - l0, s + l0, s - l1, s + l1]
    else:
        r = l0
        th = [s - l1, s + l1]
    base = s - r + sigma * (1 - a) / 2
    bad = []
    for i in range(len(th)):
        for j in range(len(th)):
            if i != j and ((th[i] - th[j]) / sigma).denominator == 1:
                bad.append(("diff", i, j))
        if ((th[i] - base) / sigma).denominator == 1:
            bad.append(("frozen", i))
    return th, base, bad


def main():
    print("== S0 table (a = 1..6 per solution; Sol 2 listed for completeness, deg L = 2 -> OUT)", flush=True)
    for sol in ("Sol3", "Sol1", "Sol2"):
        for a in range(1, 7):
            t = t_of(sol, a)
            pr = params(sol, t)
            assert pr["a"] == a, (sol, a, pr["a"])
            D = pr["degP"] + a * pr["degL"]
            pconv = F(a + 1, pr["degP"] + a) if pr["degL"] == 1 else None
            ins = pr["sigma"] > 0 and pr["degL"] == 1 and pr["degP"] >= 2
            g0 = (pr["s"]) / pr["sigma"] + F(1 - a, 2)
            print(f"  {sol} a={a}: t={t} n={pr['n']} s={pr['s']} sigma={pr['sigma']} c={pr['c']} b={pr['b']} "
                  f"degP={pr['degP']} degL={pr['degL']} | ODE order {pr['degP']}, Gamma order {D}, p={pconv}, "
                  f"q_ODE={pr['c'] / (pr['degP'] - 1)}, q_Gamma={pr['c'] / D}, g(r=0 part)={g0} | IN SCOPE: {ins} | "
                  f"record Q data: {DATA.get((sol, t), 'none')}", flush=True)
    print("== record data points not in scope (non-integer a or deg L = 2)", flush=True)
    for (sol, t), src in DATA.items():
        pr = params(sol, t)
        print(f"  {sol} t={t}: a={pr['a']}, sigma={pr['sigma']}, degL={pr['degL']}: {src}", flush=True)
    print("== (H) at the proposed S3 momentum points (l0, l1) = A (1/3, 3/5), B (9/10, 1/5)", flush=True)
    for sol, a in (("Sol3", 2), ("Sol3", 1), ("Sol3", 3), ("Sol1", 1), ("Sol1", 2)):
        pr = params(sol, t_of(sol, a))
        for nm, (l0, l1) in (("A", (F(1, 3), F(3, 5))), ("B", (F(9, 10), F(1, 5)))):
            th, base, bad = check_H(sol, pr, l0, l1)
            print(f"  {sol} a={a} t={t_of(sol, a)} point {nm}: theta_j={[str(x) for x in th]}, frozen base={base} "
                  f"(lattice base + {pr['sigma']} Z); (H) {'HOLDS' if not bad else 'FAILS ' + str(bad)}", flush=True)


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] != "S0":
        raise SystemExit("usage: spec1_s0.py S0")
    import hashlib
    import os
    me = os.path.abspath(__file__)
    print(f"script {me} sha256 {hashlib.sha256(open(me, 'rb').read()).hexdigest()} (imports: stdlib only)", flush=True)
    main()
    print("S0 DONE", flush=True)
