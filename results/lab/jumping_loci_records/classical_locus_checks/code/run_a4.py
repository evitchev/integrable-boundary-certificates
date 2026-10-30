"""Check A4: the per-sheet maps (b(t), c(t)) from the committed I_3, and their consequences.

(1) The committed I_3 densities (repository code/mixed_engine.py cyl_P4, byte copy in code/vendor, sha256-checked against the
    repository file) are evaluated over the field Q(t) (sympy FracField: canonical elements, so cyl_P4's on-curve assert is an
    exact identity).  The weight-4 symbol is the part with only order-one derivatives: (dX)^4 -> x^4, (dX)^2 (dY.dY) -> x^2 r^2,
    (dY.dY)^2 -> r^4 (code/cyl_classical_hierarchy.py lines 6-13, 279).  b = [x^2 r^2]/[x^4], c = [r^4]/[x^4], compared EXACTLY
    with the reviewers' maps.  The same for the weight-4 profiles PINS (code/cyl_mo2_profiles.py, a_1 = b, a_2 = c).
(2) Consequences: the curves c = b/6, c = b^2/(6-b), c = 1; in x = beta^2; the sheet-wise solution of eq. (3) and of the points
    at infinity; the odd-q statements; enumeration over odd s <= 61; cross-check with the repository's classical records
    (P3b EXPECT_JUMPS at weights 6-10, transcribed from code/cyl_classical_hierarchy.py lines 153-154, and the pinned lists
    results/lab/classical/rank_drop_locus_w12.txt, _w14_18.txt, _w20_22.txt, parsed)."""
import os, sys, re, hashlib
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "vendor"))
from fractions import Fraction as Fr
import sympy as sp
from sympy import QQ
import jl_lib as L

REPO = os.path.abspath(os.path.join(HERE, "..", "..", "..", "..", ".."))   # repository root (code/ -> classical_locus_checks/ -> jumping_loci_records/ -> lab/ -> results/ -> root)
fails = []
def require(cond, msg):
    print(("PASS  " if cond else "FAIL  ") + msg, flush=True)
    if not cond:
        fails.append(msg)

# ---- vendor integrity (fail-closed)
for f in ("mixed_engine.py", "invariant_engine.py", "cyl_mo2_profiles.py"):
    h1 = hashlib.sha256(open(os.path.join(HERE, "vendor", f), "rb").read()).hexdigest()
    h2 = hashlib.sha256(open(os.path.join(REPO, "code", f), "rb").read()).hexdigest()
    require(h1 == h2, f"vendor/{f} is a byte copy of the repository's code/{f} (sha256 {h1[:16]})")
if fails:
    sys.exit(1)
import mixed_engine as me
import cyl_mo2_profiles as mo2

K, T = sp.field("t", QQ)
def toK(v):
    if isinstance(v, Fr):
        return K(v.numerator) / K(v.denominator)
    if isinstance(v, int):
        return K(v)
    return v
N = (T**2 - 25) / (T**2 - 1); S = -24 * T / (T**2 - 1)
CLAIM = {1: (6 * (T - 1) / (T - 9), (T - 1) / (T - 9)),
         2: (12 * (T - 1) / (5 * T - 11), 8 * (T - 1)**2 / ((T - 3) * (5 * T - 11))),
         3: (6 * (T - 3) / (T + 1), K(1))}
KEYS = (me._M((1, 1, 1, 1)), me._M((1, 1), [(1, 1)]), me._M((), [(1, 1), (1, 1)]))

print("== (1) the maps from the committed I_3 (cyl_P4 over Q(t))")
TOP = {}
for sol in (1, 2, 3):
    P = me.cyl_P4(sol, N, S)
    X4, X2R2, R4 = (toK(P.get(k, 0)) for k in KEYS)
    TOP[sol] = (X4, X2R2, R4)
    bb, cc = X2R2 / X4, R4 / X4
    require(bb == CLAIM[sol][0], f"sheet {sol}: b(t) = [x^2 r^2]/[x^4] of the committed I_3 equals the claimed {sp.factor(CLAIM[sol][0].as_expr())}")
    require(cc == CLAIM[sol][1], f"sheet {sol}: c(t) = [r^4]/[x^4] of the committed I_3 equals the claimed {sp.factor(CLAIM[sol][1].as_expr())}")
    print(f"      sheet {sol}: the unnormalized symbol [x^4 : x^2 r^2 : r^4] = [{sp.factor(X4.as_expr())} : {sp.factor(X2R2.as_expr())} : {sp.factor(R4.as_expr())}]", flush=True)

print("== (1b) the weight-4 profiles PINS (cyl_mo2_profiles): a_1 = b, a_2 = c")
for sol in (1, 2, 3):
    prof = mo2.PINS[sol](4)
    try:
        a1, a2 = toK(prof[1](T)), toK(prof[2](T))
        route = "over Q(t)"
        ok = (a1 == CLAIM[sol][0]) and (a2 == CLAIM[sol][1])
    except Exception as e:      # the Fraction-based lambdas may refuse field arguments: then 41 rational points (degree bound below)
        route = f"at 41 rational points ({type(e).__name__} over Q(t))"
        pts = [Fr(p, q) for p, q in ((2, 1), (4, 1), (5, 1), (6, 1), (7, 1), (8, 1), (10, 1), (12, 1), (-2, 1), (-4, 1), (-5, 1),
               (-6, 1), (-7, 1), (-8, 1), (1, 2), (3, 2), (5, 2), (7, 2), (-3, 2), (-5, 2), (1, 3), (2, 3), (4, 3), (5, 3), (7, 3),
               (-2, 3), (-4, 3), (1, 4), (3, 4), (5, 4), (7, 4), (-3, 4), (1, 5), (2, 5), (3, 5), (4, 5), (6, 5), (7, 5), (8, 5),
               (9, 7), (13, 11))]
        ok = True
        for p in pts:
            num = sp.Rational(p.numerator, p.denominator)
            for k_, claim in ((1, CLAIM[sol][0]), (2, CLAIM[sol][1])):
                cv = claim.as_expr().subs(sp.Symbol('t'), num)
                ok &= sp.Rational(prof[k_](p).numerator, prof[k_](p).denominator) == cv
        # each side: numerator/denominator of degree <= 2; the difference cleared has degree <= 4 < 41
    require(ok, f"sheet {sol}: PINS[{sol}](4) = [1, b, c] with the claimed b, c ({route})")

print("== (2) consequences")
x = sp.Symbol('x'); ts = sp.Symbol('t'); bs, cs = sp.symbols('b c')
def claim_expr(sol):
    return tuple(sp.together(v.as_expr()) for v in CLAIM[sol])
in_x = {}
for sol in (1, 2, 3):
    bb, cc = claim_expr(sol)
    bx = sp.factor(sp.cancel(bb.subs(ts, (1 + x) / (1 - x)))); cx = sp.factor(sp.cancel(cc.subs(ts, (1 + x) / (1 - x))))
    in_x[sol] = (bx, cx)
    print(f"      sheet {sol} in x = beta^2: b = {bx}, c = {cx}", flush=True)
require(sp.cancel(in_x[1][1] - in_x[1][0] / 6) == 0, "sheet 1: c = b/6 identically")
require(sp.cancel(in_x[2][1] - in_x[2][0]**2 / (6 - in_x[2][0])) == 0, "sheet 2: c = b^2/(6 - b) identically")
require(sp.expand(in_x[3][0] - (12 * x - 6)) == 0 and in_x[3][1] == 1, "sheet 3: b = 12 beta^2 - 6, c = 1")

# eq. (3) on each sheet, symbolically in (s, j, j')
sS, jj, jp = sp.symbols('s j jp', positive=True)
b3 = 6 * jj / (3 * jj - sS); c3 = jj * (3 * jp - sS) / (jp * (3 * jj - sS))
require(sp.solve(sp.Eq(c3, b3 / 6), sS) == [2 * jp], "sheet 1 (c = b/6) on eq. (3): s = 2j' (even): no affine classical point at odd spin")
require(sp.solve(sp.Eq(sp.Rational(1, 6), (3 * jp - sS) / (6 * jp)), sS) == [2 * jp],
        "sheet 1 at infinity (its normalization pole t = 9: c/b = 1/6) against C' = (3j'-s)/(6j'): s = 2j' (even): none at odd spin")
sol2 = sp.solve(sp.Eq(c3 * (6 - b3), b3**2), sS)
require(set(sol2) == {2 * jj + 3 * jp}, f"sheet 2 (c = b^2/(6-b)) on eq. (3): s = 2j + 3j' (solutions {sol2})")
require(sp.solve(sp.Eq(c3, 1), jp) == [jj], "sheet 3 (c = 1) on eq. (3): j' = j")
x3 = sp.solve(sp.Eq(12 * x - 6, b3.subs(jp, jj)), x)
require(len(x3) == 1 and sp.cancel(x3[0] - (sS - 4 * jj) / (2 * (sS - 3 * jj))) == 0, f"sheet 3: beta^2 = (s - 4j)/(2(s - 3j)) ({x3})")
x2 = sp.solve(sp.Eq(in_x[2][0], b3.subs(sS, 2 * jj + 3 * jp)), x)
require(len(x2) == 1 and sp.cancel(x2[0] - jj / (2 * (jj + jp))) == 0, f"sheet 2: beta^2 = j/(2(j + j')) (the hand derivation) ({x2})")

# enumeration over odd s <= 61
SMAX = 61
odd_q_sheet3 = []; even_q_sheet3 = []; odd_q_sheet2 = []; even_q_sheet2 = []; bad = []
for s in range(3, SMAX + 1, 2):
    for j, tv in L.sheet3_points(s).items():
        xv = L.x_of(tv)
        if xv != Fr(s - 4 * j, 2 * (s - 3 * j)):
            bad.append(("sheet3 x", s, j))
        b_, c_ = L.lattice_affine(s)[(j, j)]
        bt, ct = (sp.Rational(v.as_expr().subs(sp.Symbol('t'), sp.Rational(tv.numerator, tv.denominator))) for v in CLAIM[3])
        if (bt, ct) != (sp.Rational(b_.numerator, b_.denominator), sp.Rational(c_.numerator, c_.denominator)):
            bad.append(("sheet3 map", s, j))
        if xv.numerator % 2 == 0 or xv.denominator % 2 == 1:
            bad.append(("sheet3 parity", s, j, xv))
        if xv.numerator == 1:
            (even_q_sheet3 if xv.denominator % 2 == 0 else odd_q_sheet3).append((s, xv.denominator, j, tv))
    for (j, jp_), tv in L.sheet2_points(s).items():
        xv = L.x_of(tv)
        if xv != Fr(j, 2 * (j + jp_)):
            bad.append(("sheet2 x", s, j, jp_))
        if j != 3 * jp_:
            b_, c_ = L.lattice_affine(s)[(j, jp_)]
            bt, ct = (sp.Rational(v.as_expr().subs(sp.Symbol('t'), sp.Rational(tv.numerator, tv.denominator))) for v in CLAIM[2])
            if (bt, ct) != (sp.Rational(b_.numerator, b_.denominator), sp.Rational(c_.numerator, c_.denominator)):
                bad.append(("sheet2 map", s, j, jp_))
        else:
            if L.lattice_infinity(s).get((j, jp_)) != Fr(-1):
                bad.append(("sheet2 infinity", s, j, jp_))
        if xv.numerator == 1:
            (even_q_sheet2 if xv.denominator % 2 == 0 else odd_q_sheet2).append((s, xv.denominator, j, jp_, tv))
require(not bad, f"odd s <= {SMAX}: every closed-form point lies on its sheet's curve at the stated t (eq. (3) or the point at infinity); "
        f"sheet-3 beta^2 always odd/even in lowest terms ({len(bad)} exceptions{': ' + str(bad[:5]) if bad else ''})")
require(not odd_q_sheet3, f"sheet 3: beta^2 = 1/q with q odd never occurs classically (odd s <= {SMAX}); q even occurs {len(even_q_sheet3)} times")
print("      sheet 3, beta^2 = 1/q (q even), s <= 61: " + ", ".join(f"s={s} q={q} t={tv}" for s, q, j, tv in even_q_sheet3), flush=True)
print("      sheet 2, beta^2 = 1/q (q odd),  s <= 61: " + ", ".join(f"s={s} q={q} (j,j')=({j},{jp_}) t={tv}" for s, q, j, jp_, tv in odd_q_sheet2), flush=True)
print("      sheet 2, beta^2 = 1/q (q even), s <= 61: " + ", ".join(f"s={s} q={q} (j,j')=({j},{jp_}) t={tv}" for s, q, j, jp_, tv in even_q_sheet2), flush=True)
pattern = all(s % (3 * q - 2) == 0 and (s // (3 * q - 2)) % 2 == 1 and j == 2 * (s // (3 * q - 2)) and jp_ == (q - 2) * (s // (3 * q - 2))
              for s, q, j, jp_, tv in odd_q_sheet2)
expect_pat = {(m * (3 * q - 2), q) for q in range(3, 41, 2) for m in range(1, 41, 2) if m * (3 * q - 2) <= SMAX}
require(pattern and {(s, q) for s, q, *_ in odd_q_sheet2} == expect_pat,
        "sheet 2: beta^2 = 1/q with q odd occurs exactly at s = m(3q - 2), m odd, with (j, j') = (2m, m(q-2)) (s <= 61)")
lowest = sorted((s, q, tv) for s, q, j, jp_, tv in odd_q_sheet2 if j == 2)
require([(s, tv) for s, q, tv in lowest if q in (3, 5, 7)] == [(7, Fr(2)), (13, Fr(3, 2)), (19, Fr(4, 3))],
        "sheet 2: the lowest odd-q points (j, j') = (2, q-2): t = 2, 3/2, 4/3 at spins 7, 13, 19")

print("== (3) the repository's classical lists against the closed forms")
EXPECT_JUMPS = {3: {6: {Fr(5, 3)}, 8: {Fr(11, 5)}, 10: {Fr(7, 5), Fr(17, 7)}},     # cyl_classical_hierarchy.py lines 153-154 (P3b,
                2: {6: {Fr(5, 3)}, 8: {Fr(2)}, 10: set()}, 1: {6: set(), 8: set(), 10: set()}}   # the exhaustive affine enumeration)
def parse_family(path):
    out = {1: {}, 2: {}, 3: {}}; Wc = None
    for line in open(path):
        m = re.match(r'weight (\d+):', line)
        if m: Wc = int(m.group(1))
        m = re.match(r'   Solution (\d): (?:non-decoupling )?classical jumps at t = \[(.*?)\]', line)
        if m and Wc:
            out[int(m.group(1))][Wc] = {Fr(v) for v in re.findall(r'-?\d+(?:/\d+)?', m.group(2))}
    return out
pinned = {1: {}, 2: {}, 3: {}}
for fn in ("rank_drop_locus_w12.txt", "rank_drop_locus_w14_18.txt", "rank_drop_locus_w20_22.txt"):
    pf = parse_family(os.path.join(REPO, "results", "lab", "classical", fn))
    for sol in (1, 2, 3):
        pinned[sol].update(pf[sol])
require(all(set(pinned[sol]) == {12, 14, 16, 18, 20, 22} for sol in (1, 2, 3)), "the pinned lists cover weights 12-22 for the three sheets")
for Wt in range(6, 24, 2):
    s = Wt - 1
    src = EXPECT_JUMPS if Wt <= 10 else pinned
    mine3 = set(L.sheet3_points(s).values()); mine2_aff = {tv for (j, jp_), tv in L.sheet2_points(s).items() if j != 3 * jp_}
    mine2_inf = {tv for (j, jp_), tv in L.sheet2_points(s).items() if j == 3 * jp_}
    rec = {sol: src[sol][Wt] for sol in (1, 2, 3)}
    require(rec[3] == mine3 and rec[2] == mine2_aff and rec[1] == set(),
            f"weight {Wt} (spin {s}): sheet 3 {sorted(mine3)}, sheet 2 affine {sorted(mine2_aff)}, sheet 1 none: equal to "
            f"{'P3b EXPECT_JUMPS' if Wt <= 10 else 'the pinned list'}; sheet 2 at infinity: {sorted(mine2_inf) if mine2_inf else 'none'}")
for Wt, tv in ((8, Fr(2)), (14, Fr(3, 2)), (20, Fr(4, 3))):
    src = EXPECT_JUMPS if Wt <= 10 else pinned
    require(tv in src[2][Wt], f"the odd-q point t = {tv} is in the repository's sheet-2 list at weight {Wt}")

print()
if fails:
    print(f"A4: {len(fails)} FAILURE(S)"); sys.exit(1)
print("A4 DONE: all checks pass")
