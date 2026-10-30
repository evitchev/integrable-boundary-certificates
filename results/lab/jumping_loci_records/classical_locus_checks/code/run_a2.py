"""Check A2: eq. (3) needs 3j != s; the projective stratum; sheet 2 at t = 11/5.

(1) Admissible lattice indices with 3j = s (s = 3..29): eq. (3) is undefined there.
(2) The line at infinity al = 0, be = 1: H = x^2 r^2 + C' r^4, C' = c/b.  For W = 4..30: all maximal minors of the conservation
    matrix as polynomials in C' (determinants computed, not assumed), their gcd, its roots; compared with the set
    {(3j'-s)/(6j') : 3j = s, j + j' <= (s-1)/2}; the nullity at each root (exact) and at non-roots.
(3) The decoupling line be = 0 (rank over Q(al, ga)) and the point r^4.
(4) Which sheets meet the line at infinity at a fiber (al(t) = 0 of the homogeneous maps), and where.
(5) Sheet 2 at t = 11/5 from the committed I_3 (cyl_P4 at (N, s) = (-21/4, -55/4), byte copy, Fractions): its weight-4 symbol
    and its nullity at W = 4..30, against the certificate's P3e (code/cyl_classical_hierarchy.py lines 265-283)."""
import os, sys, hashlib
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "vendor"))
from fractions import Fraction as Fr
import sympy as sp
import jl_lib as L
from jl_lib import x, r, al, be, ga, t

REPO = os.path.abspath(os.path.join(HERE, "..", "..", "..", "..", ".."))   # repository root (code/ -> classical_locus_checks/ -> jumping_loci_records/ -> lab/ -> results/ -> root)
fails = []
def require(cond, msg):
    print(("PASS  " if cond else "FAIL  ") + msg, flush=True)
    if not cond:
        fails.append(msg)

def rat(v):
    return sp.Rational(v.numerator, v.denominator) if isinstance(v, Fr) else sp.Rational(v)

def nul_fr(W, A, B, C):
    return L.nullity_at(W, rat(A), rat(B), rat(C))

print("== (1) admissible indices with 3j = s")
und = {s: [(j, jp) for j in range(1, s) for jp in range(1, s) if j + jp <= (s - 1) // 2 and 3 * j == s] for s in range(3, 31, 2)}
for s, lst in und.items():
    if lst:
        print(f"      s={s}: (j, j') = {lst}: 3j - s = 0, eq. (3) undefined", flush=True)
require(und[9] == [(3, 1)], "s = 9, j = 3, j' = 1 is admissible (j + j' = 4 = (s-1)/2) and makes eq. (3) undefined")

print("== (2) the line at infinity: H = x^2 r^2 + C' r^4")
cp = sp.Symbol('Cp')
inf_sets = {}
for W in range(4, 32, 2):
    k = W // 2; s = W - 1
    M = L.matrix(W, L.quartic(0, 1, cp))
    minors = [sp.expand(M.extract(list(range(k)), [j for j in range(k + 1) if j != m]).det(method='berkowitz')) for m in range(k + 1)]
    g = sp.Integer(0)
    for mn in minors:
        g = sp.gcd(g, mn)
    g = sp.Poly(g, cp)
    roots = []
    nonlinear = []
    if g.degree() > 0:
        for fac, _ in sp.factor_list(g.as_expr())[1]:
            pf = sp.Poly(fac, cp)
            if pf.degree() == 1:
                cf = pf.all_coeffs(); roots.append(Fr(int(sp.Rational(-cf[1] / cf[0]).p), int(sp.Rational(-cf[1] / cf[0]).q)))
            elif pf.degree() > 1:
                nonlinear.append(fac)
    expect = set(L.lattice_infinity(s).values())
    ok = not nonlinear and set(roots) == expect
    nul_ok = all(nul_fr(W, 0, 1, rv) == 2 for rv in roots)
    ctrl = [Fr(1, 6), Fr(0), Fr(1), Fr(-1, 2), Fr(5, 7)]
    ctrl_ok = all(nul_fr(W, 0, 1, cv) == 1 for cv in ctrl if cv not in expect)
    inf_sets[W] = sorted(roots)
    require(ok and nul_ok and ctrl_ok, f"W={W} (spin {s}): gcd of the {k + 1} maximal minors has degree {g.degree()}; drop values C' = "
            f"{[str(v) for v in sorted(roots)]} = the 3j = s lattice {[str(v) for v in sorted(expect)]}; nullity 2 at each, 1 at the controls")
require(inf_sets[10] == [Fr(-1)] and all(inf_sets[W] == [] for W in (6, 8, 12, 14)),
        "P3e's statement for weights 6-14 (only c/b = -1, at weight 10) holds")

print("== (3) the decoupling line be = 0")
for W in (6, 10, 16, 22, 30):
    k = W // 2
    M0 = L.matrix(W, L.quartic(al, 0, ga))
    rk = M0.rank()
    require(M0.cols - rk == 2, f"W={W}: at be = 0 the matrix has rank {rk} over Q(al, ga): nullity exactly 2 (F_xr = 0: A x^W + B r^W)")

print("== (4) the sheets against the line at infinity")
for sheet in (1, 2, 3):
    A, B, C = L.homogeneous(sheet)
    zs = [z for z in sp.solve(A, t)]
    info = []
    for z in zs:
        Bz, Cz = sp.simplify(B.subs(t, z)), sp.simplify(C.subs(t, z))
        info.append((z, sp.Rational(Cz / Bz) if Bz != 0 else None))
    print(f"      sheet {sheet}: homogeneous [al : be : ga] = [{sp.factor(A)} : {sp.factor(B)} : {sp.factor(C)}]; al = 0 at t = {zs}: (t, C') = {info}", flush=True)
    for z, cpv in info:
        if z in (1, -1):
            continue
        if cpv is None:     # be(z) = 0 as well: the point [0 : 0 : 1] = r^4, on the decoupling line
            require(sheet == 2 and z == 3, f"sheet {sheet} at t = {z}: the point r^4 = [0 : 0 : 1] on the decoupling line be = 0 (A6)")
            continue
        hits = [W for W in range(4, 32, 2) if cpv is not None and Fr(int(cpv.p), int(cpv.q)) in set(inf_sets[W])]
        if sheet == 2:
            require(z == sp.Rational(11, 5) and cpv == -1 and hits == [10, 28],
                    f"sheet 2 meets the line at infinity at t = {z} (a fiber; beta^2 = {L.x_of(Fr(11, 5))}), C' = {cpv}: a drop exactly at W = {hits} (spins {[W - 1 for W in hits]})")
        else:
            require(hits == [], f"sheet {sheet} meets the line at infinity at t = {z}, C' = {cpv}: no drop at any W = 4..30")

print("== (5) sheet 2 at t = 11/5 from the committed I_3")
for f in ("mixed_engine.py", "invariant_engine.py"):
    h1 = hashlib.sha256(open(os.path.join(HERE, "vendor", f), "rb").read()).hexdigest()
    h2 = hashlib.sha256(open(os.path.join(REPO, "code", f), "rb").read()).hexdigest()
    require(h1 == h2, f"vendor/{f} is a byte copy of the repository's code/{f}")
import mixed_engine as me
tv = Fr(11, 5); Nv, sv = L.N_of(tv), L.s_of(tv)
require((Nv, sv) == (Fr(-21, 4), Fr(-55, 4)), f"t = 11/5 is the fiber (N, s) = ({Nv}, {sv})")
P = me.cyl_P4(2, Nv, sv)
X4, X2R2, R4 = (P.get(kk, Fr(0)) for kk in (me._M((1, 1, 1, 1)), me._M((1, 1), [(1, 1)]), me._M((), [(1, 1), (1, 1)])))
require((X4, X2R2, R4) == (Fr(0), Fr(-1), Fr(1)), f"the committed I_3 of sheet 2 at t = 11/5 has weight-4 symbol [x^4 : x^2 r^2 : r^4] = [{X4} : {X2R2} : {R4}] "
        f"= -r^2 (x^2 - r^2): the x^4 coefficient vanishes (the affine chart fails), C' = c/b = {R4 / X2R2}")
nul = {W: nul_fr(W, X4, X2R2, R4) for W in range(4, 32, 2)}
require(all(nul[W] == (2 if (W - 1) % 9 == 0 else 1) for W in nul),
        f"its nullity at W = 4..30: {nul}: 2 exactly at W = 10, 28 (spins 9, 27 = 9(2h+1)), 1 elsewhere -- P3e's 'nullity 2 at weights 10 and 28, 1 at 6, 8, 12' extended")
print()
if fails:
    print(f"A2: {len(fails)} FAILURE(S)"); sys.exit(1)
print("A2 DONE: all checks pass")
