"""Exact-elimination certificate (2026-09-04): at the sampled pole positions, a
momentum-dependent fourth regular singular point is excluded EXACTLY under
the weaker (proportionality) reading of rationality, for the m = -1 Heun
family and for the sliding-pole family.

Setting.  Heun family (Solution 3, m = -1):
    Psi'' = [kappa^2/(1+u) + A u/(1+u) + B u/(1+u)^2 + D u/(u+v) + E u/(u+v)^2] Psi;
sliding-pole family: leading symbol kappa^2/(u+v), the same four terms.
At a fixed rational pole position v the WKB coefficient S_s is
R0 + T R1 with one transcendental T (asin-type at v > 1); the certified
VEVs are rational in t, so with the correspondence S_s = alpha_s I_s + beta_s
and alpha_s allowed to contain T (the physically natural reading: BLZ-type
normalizations are transcendental in the family parameter) the rational
part and the T-part of S_s must EACH be zero or proportional to I_s
(1 and T being linearly independent over Q(t) for a non-constant algebraic
v(t)).  At top degree in the momenta that is: for every spin s and each
part, the top form lies in {0} u C * profile_s, where profile_s is the
certified palindromic Solution-3 profile (u_4 = w free; u_6 = 15w/(12-w),
u_8 = 14w/(9-w), m_8 = 42w(2w+3)/((9-w)(w+12)) from the certified closed
forms).  A momentum-dependent fourth point means some top coefficient of
D or E is nonzero: scaled to 1 in four branches; w not in {+-2, 6}, i.e.
t not in {5, 2, infinity} (the two single-form degeneracies and the
N = 1 anchor, where the fourth point must disappear), imposed by a
Rabinowitsch variable.  A Groebner basis {1} in every branch = no such
operator, exactly.

Inputs (PINNED): the general S_1..S_5 (S_1..S_7 for the sliding pole) at
each position as polynomials in (A, B, D, E) over Q, split into the
rational part and the T-coefficient -- results/lab/heun/general_S/*.txt,
produced by lab/heun_general_S.py and lab/slide_general_S.py from the lab
integrators (Euler/2F1 continuation, checked against numerical
quadrature); the spin-1 D-term of the Heun family is re-derived here
by hand (P0) as a check of the pinned files' bookkeeping.

Checks:
  P0  the pinned Heun S_1 at each position has the form
      -A + B - D + E/(2v) - 1/8 - T[2D(v-1) - E]/(2v);
  P1  Heun family, v = 3, 7/3, 5/2, 9, 1/2, 2: Groebner basis {1} in all
      four branches (spins 1-5; at v = 2 the exclusion of w = 6 is needed,
      elsewhere w != +-2 already suffices);
  P2  sliding pole, v = 3, spins 1-7: basis {1} in all four branches;
  P3  negative control: with the fourth point removed (D = E = 0 forced,
      u_4 = 6 allowed) the same machinery returns a NON-trivial basis --
      the two-term case has its solutions at the N = 1 degeneracy -- so a
      basis {1} is not what the code returns regardless.
Fail-closed: any failure exits 1 without the CERTIFIED marker.
Runtime ~2-3 min.
"""
import sys, os, time
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
GS = os.path.join(ROOT, "results", "lab", "heun", "general_S")
fails = []
def require(cond, msg):
    print(("PASS  " if cond else "FAIL  ") + msg, flush=True)
    if not cond:
        fails.append(msg)

A, B, D, E, T = sp.symbols('A B D E T')
x, y, w, rb = sp.symbols('x y w rb')
a1, a2, b1, b2, d1, d2, e1, e2 = cs = sp.symbols('a1 a2 b1 b2 d1 d2 e1 e2')
SUB = {A: a1*x + a2*y, B: b1*x + b2*y, D: d1*x + d2*y, E: e1*x + e2*y}

def load_S(fname):
    S = {}
    with open(os.path.join(GS, fname)) as fh:
        for line in fh:
            name, expr = line.split(' = ', 1); S[name] = eval(expr, sp.__dict__)
    return S

def top(expr, deg):
    P_ = sp.Poly(expr, A, B, D, E)
    return sum(cf * A**i*B**j*D**k*E**l for (i, j, k, l), cf in P_.terms() if i + j + k + l == deg)

def prof_eqs(form, deg, palindromic=True):
    P_ = sp.Poly(sp.expand(form), x, y); c = lambda i, j: P_.coeff_monomial(x**i * y**j)
    if deg == 1: return [c(1, 0) - c(0, 1)] if palindromic else []
    if deg == 2: return [c(2, 0) - c(0, 2), c(1, 1) - w * c(2, 0)] if palindromic else [c(1, 1) - w * c(2, 0)]
    if deg == 3: return [c(3, 0) - c(0, 3), c(2, 1) - c(1, 2), (12 - w) * c(2, 1) - 15 * w * c(3, 0)]
    if deg == 4: return [c(4, 0) - c(0, 4), c(3, 1) - c(1, 3), (9 - w) * c(3, 1) - 14 * w * c(4, 0), (9 - w) * (w + 12) * c(2, 2) - 42 * w * (2 * w + 3) * c(4, 0)]

def basis(S, kmax, branch, exclude6=False, palindromic=True, force_no_point=False):
    eqs = []
    for k in range(1, kmax + 1, 2):
        deg = (k + 1) // 2
        for tag in ('T0', 'T1'):
            eqs += prof_eqs(top(S[f'S{k}_{tag}'], deg).subs(SUB), deg, palindromic)
    scale = {'d1': d1, 'd2': d2, 'e1': e1, 'e2': e2}[branch]
    excl = (w**2 - 4) * ((w - 6) if exclude6 else 1)
    if force_no_point:
        eqs = [sp.expand(e_.subs({d1: 0, d2: 0, e1: 0, e2: 0, a1: 1})) for e_ in eqs] + [1 - rb * excl]
        gens = [a2, b1, b2, w, rb]
    else:
        eqs = [sp.expand(e_.subs(scale, 1)) for e_ in eqs] + [1 - rb * excl]
        gens = [c_ for c_ in cs if c_ != scale] + [w, rb]
    eqs = [sp.numer(sp.together(e_)) for e_ in eqs if e_ != 0]
    return list(sp.groebner(eqs, *gens, order='grevlex'))

if __name__ == "__main__":
    t0 = time.time()
    heun = {"3": ("heun_S_v3_k5.txt", sp.Integer(3)), "7/3": ("heun_S_v7_3_k5.txt", sp.Rational(7, 3)), "5/2": ("heun_S_v5_2_k5.txt", sp.Rational(5, 2)),
            "9": ("heun_S_v9_k7.txt", sp.Integer(9)), "1/2": ("heun_S_v1_2_k5.txt", sp.Rational(1, 2)), "2": ("heun_S_v2_k5.txt", sp.Integer(2))}
    # P0: the pinned spin-1 bookkeeping against the hand derivation
    ok0 = True
    for label, (fname, vv) in heun.items():
        S = load_S(fname)
        s1 = sp.expand(S['S1_T0'] + T * S['S1_T1'])
        hand = sp.expand(-A + B - D + E/(2*vv) - sp.Rational(1, 8) - T * (2*D*(vv - 1) - E)/(2*vv))
        ok0 &= sp.expand(s1 - hand) == 0
    require(ok0, "P0: the pinned Heun S_1 at all six positions equals -A + B - D + E/(2v) - 1/8 - T[2D(v-1) - E]/(2v) (hand derivation)")
    # P1: Heun family
    for label, (fname, vv) in heun.items():
        S = load_S(fname); res = {}
        for br in ('d1', 'd2', 'e1', 'e2'):
            res[br] = basis(S, 5, br, exclude6=True)
        ok = all(res[br] == [1] for br in res)
        require(ok, f"P1: Heun family at v = {label}: Groebner basis {{1}} in all four scaling branches (spins 1-5, weak reading, u_4 not in {{+-2, 6}})"
                + ("" if ok else f" -- {[(br, str(res[br])[:60]) for br in res if res[br] != [1]]}"))
    # P2: sliding pole, spins 1-7
    S = load_S("slide_S_v3_k7.txt"); res = {br: basis(S, 7, br, exclude6=True) for br in ('d1', 'd2', 'e1', 'e2')}
    ok = all(res[br] == [1] for br in res)
    require(ok, "P2: sliding pole at v = 3: Groebner basis {1} in all four branches (spins 1-7, weak reading, u_4 not in {+-2, 6})"
            + ("" if ok else f" -- {[(br, str(res[br])[:60]) for br in res if res[br] != [1]]}"))
    # P3: negative controls
    S = load_S("heun_S_v3_k5.txt")
    two = basis(S, 5, 'd1', force_no_point=True)
    require(two != [1] and all(str(g_) for g_ in two), "P3: negative control -- with the fourth point removed the system is the two-term case and has solutions (the u_4 = +-2, 6 degeneracies)")
    print(f"\n{time.time()-t0:.0f}s")
    if fails:
        print(f"EXACT ELIMINATION: {len(fails)} FAILURE(S)"); sys.exit(1)
    print("EXACT ELIMINATION CERTIFIED (weak reading of rationality: no momentum-dependent fourth point for the m = -1 Heun family at "
          "v = 3, 7/3, 5/2, 9, 1/2, 2 nor for the sliding pole at v = 3, from the top-degree conditions of spins 1-5 / 1-7, exactly)")
