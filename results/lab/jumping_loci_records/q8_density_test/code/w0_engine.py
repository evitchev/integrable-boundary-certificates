"""q = 8 test, w0: ENGINE CONTROLS (C1, C7, C8, C9). (a) kappa = 1 reproduces WUW1; (b) kappa = 3/2 reproduces
WUW2 (I_3 gate, W^(5), P2's unique bilinear, dim K_10 = 2 by the Euler + fmpq_mat route on C_10, the eigenvalue gate, P4b); (c) Euler == Ech
and rank(E) = p2(w) - p2(w-1); C_w checks. Exit 0 iff every control passes."""
import os, sys, time, random, hashlib
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import q8_lib as L
from q8_lib import K, A, B, say
import sympy as sp
say('w0_engine.py sha256', L.self_hash(__file__))
T0 = time.time(); FAIL = []
def control(name, ok):
    ok = bool(ok); say(f'CONTROL {name}: {ok}  ({time.time()-T0:.0f}s)')
    if not ok: FAIL.append(name)
rng = random.Random(20261028)
def rand_poly(w, n=6):
    ms = K.monomials(w); return {m: F(rng.randint(-5, 5) or 1, rng.randint(1, 4)) for m in rng.sample(ms, min(n, len(ms)))}
# ---------------- (c) complement basis (C8), Euler rank (C9), Euler == Ech
for w in (6, 10, 14):
    n0, nC, new, tot, p2w, p2w1 = L.check_complement(w)
    control(f'C8 C_{w}: |dF_{w-1}| = {n0} (p2 = {p2w1}), |C_w| = {nC}, all independent mod dF: {new == nC}, dF + C_w = {tot} (p2({w}) = {p2w})', n0 == p2w1 and new == nC == p2w - p2w1 and tot == p2w)
for w in (6, 10):
    vecs = [L.euler_vec({m: F(1)}) for m in K.monomials(w)]
    rk = L.rank_Q(vecs); control(f'C9 rank(E_V, E_Z) on F_{w} = {rk} == p2(w) - p2(w-1) = {L.p2(w) - L.p2(w-1)}', rk == L.p2(w) - L.p2(w - 1))
for w in (10, 13):
    e = K.d_ech(w); ok = True; nz = 0
    for _ in range(20):
        P = rand_poly(w); Q = K.dpoly(rand_poly(w - 1)); S = K.add(P, Q)
        for X in (P, Q, S):
            ex_e = not e.reduce(dict(X)); ex_u = L.is_exact(X); ok &= (ex_e == ex_u); nz += (not ex_u)
    control(f'C9 Euler == Ech at weight {w} (20 random P, dQ, P + dQ; non-exact seen {nz}/60)', ok and 0 < nz < 60)
# ---------------- (a) kappa = 1: WUW1
K.set_kappa(1)
c4 = K.reps_mod_d(K.commutant(4, [{A: F(1, 2), B: F(1, 2)}, {A: F(2), B: F(0)}, {A: F(0), B: F(2)}]), 4)
WU1 = K.current(2, A)
control('(a) kappa = 1: I_3 dim 1; W_U = (4/3)(dU)^3 + 2 dU d^2U + (1/3) d^3U', len(c4) == 1 and WU1 == {((A, 0),) * 3: F(4, 3), ((A, 0), (A, 1)): F(2), ((A, 2),): F(1, 3)})
fA = [{m: F(1)} for m in K.monomials(3, (A,))]; fB = [{m: F(1)} for m in K.monomials(3, (B,))]
prods = [K.mul(f, g) for f in fA for g in fB]
vecs = [L.euler_vec(K.ope(c4[0], p_, poles={1}).get(1, {})) for p_ in prods]
rel, rk = L.nullspace_Q(vecs)
wv = [WU1.get(next(iter(f)), F(0)) for f in fA]; target = [wv[i] * wv[j] for i in range(3) for j in range(3)]
k0 = next(k for k in range(9) if target[k] != 0)
control('(a) kappa = 1: the bilinear space in K_6 + d is 1-dim = W_U (x) W_W, i.e. (1, 3/2, 1/4)^2 (Euler route)', len(rel) == 1 and all(rel[0][k] * target[k0] == target[k] * rel[0][k0] for k in range(9)))
control('C7 eigen of total derivatives = 0 (kappa = 1, 5 random weight-5 fields)', all(L.eigen_raw(K.dpoly(rand_poly(5))) == 0 for _ in range(5)))
# ---------------- (b) kappa = 3/2: WUW2
K.set_kappa(F(3, 2))
c4 = K.reps_mod_d(K.commutant(4, [{A: F(1, 3), B: F(1, 3)}, {A: F(2), B: F(0)}, {A: F(0), B: F(2)}]), 4); I3 = c4[0]
cert3 = L.read_cyl34('3', '7/5')
control('(b) kappa = 3/2: I_3 dim 1 and eigen ~ cyl(3,4) at 7/5 (record vev_w4_points.json)', len(c4) == 1 and L.proportional(L.internal(L.eigen_raw(I3), F(3, 2)), cert3, (L.PX, L.PI)))
W5 = K.current(2, A); W5Z = K.current(2, B)
W5rec = {((1, 0), (1, 0), (1, 0), (1, 0), (1, 0)): F(4, 15), ((1, 0), (1, 0), (1, 0), (1, 1)): F(4, 3), ((1, 0), (1, 0), (1, 2)): F(2, 3), ((1, 0), (1, 1), (1, 1)): F(1), ((1, 0), (1, 3)): F(1, 6), ((1, 1), (1, 2)): F(1, 3), ((1, 4),): F(1, 60)}
control('(b) W^(5) equals the record (run_wuw2.log line 4)', W5 == W5rec)
eW = L.eigen_raw(W5); P = L.PA
control('(b) P1 of WUW2: eigen(W^(5)) = (1/15) P (P^2 - 1)(4P^2 - 1) = binom(2P + 2, 5) exactly', sp.expand(eW - sp.Rational(1, 15) * P * (P**2 - 1) * (4 * P**2 - 1)) == 0)
fA = [{m: F(1)} for m in K.monomials(5, (A,))]; fB = [{m: F(1)} for m in K.monomials(5, (B,))]
prods = [K.mul(f, g) for f in fA for g in fB]
vecs = [L.euler_vec(K.ope(I3, p_, poles={1}).get(1, {})) for p_ in prods]
rel, rk = L.nullspace_Q(vecs)
wv = [W5.get(next(iter(f)), F(0)) for f in fA]; target = [wv[i] * wv[j] for i in range(7) for j in range(7)]
k0 = next(k for k in range(49) if target[k] != 0)
control('(b) P2 of WUW2: the 49-bilinear space in K_10 + d is 1-dim = W^(5) (x) W^(5) (Euler route; rank %d)' % rk, len(rel) == 1 and all(rel[0][k] * target[k0] == target[k] * rel[0][k0] for k in range(49)))
C10 = L.complement_basis(10); t1 = time.time()
vecs = [L.euler_vec(K.ope(I3, {m: F(1)}, poles={1}).get(1, {})) for m in C10]
say(f'   K_10 route: {len(C10)} columns, images + Euler {time.time()-t1:.0f}s')
rel, rk = L.nullspace_Q(vecs); rkp = L.rank_mod(vecs, 2147483647)
control(f'(b) dim K_10 = {len(rel)} == 2 (record/VIR5d) by the Euler + fmpq_mat route on C_10; rank {rk} over Q, {rkp} mod 2^31-1', len(rel) == 2 and rk == rkp)
reps = [{C10[j]: x[j] for j in range(len(C10)) if x[j] != 0} for x in rel]
G10 = K.mul(W5, W5Z)
control('(b) P3 of WUW2: Euler image of Res_1 I_3 x G_10 vanishes (G_10 in K_10 + d)', L.is_exact(K.ope(I3, G10, poles={1}).get(1, {})))
eG = L.internal(L.eigen_raw(G10), F(3, 2)); p2_ = sp.Rational(1, 3)
E2 = sp.expand(sp.prod([(L.PX + L.PI - j * sp.sqrt(p2_)) * (L.PX - L.PI - j * sp.sqrt(p2_)) for j in range(-2, 3)]))
control('(b) eigen(G_10) = -(81/4800) E_2 exactly (record)', sp.expand(eG + sp.Rational(81, 4800) * E2) == 0)
vp = os.path.join(L.REPO, 'results/lab/vev/vev_sol3_w10.json')
control('(b) record file vev_sol3_w10.json has WUW2 inputs/ORIGIN_SHA256 hash 51eecbbc...', hashlib.sha256(open(vp, 'rb').read()).hexdigest().startswith('51eecbbcfb14e32c'))
cert = L.read_vev_coeffs(vp, sp.Rational(7, 5)); eK = [L.internal(L.eigen_raw(r), F(3, 2)) for r in reps]
cs = sp.symbols('c0:2'); sol = sp.solve(sp.Poly(sp.expand(cert - sum(c * e_ for c, e_ in zip(cs, eK))), L.PX, L.PI).coeffs(), cs, dict=True)
control('(b) GATE of k10_full: certified vev_sol3_w10 at 7/5 in span eigen(K_10)', bool(sol))
TT = {((A, 0), (A, 0)): F(1, 3), ((B, 0), (B, 0)): F(1, 3), ((A, 1),): F(-5, 6), ((B, 1),): F(-5, 6)}
vecs = [L.poles_vec(TT, K.dpoly({h: F(1)})) for h in K.monomials(9)]; rhs = L.poles_vec(TT, G10)
okc, r0, r1 = L.consistent_Q(vecs, rhs)
control(f'(b) P4b of WUW2 by the fmpq_mat route: no T_total-primary representative G_10 + dH (rank {r0} -> augmented {r1})', not okc and r0 == 300)
TU = {((A, 0), (A, 0)): F(1, 3), ((A, 1),): F(2, 3)}
vac = L.vir_vacuum_fields(TU, 5)
control('(b) W^(5) not in the weight-5 T_U vacuum module span (%s)' % [p for p, _ in vac], not L.in_span([f_ for _, f_ in vac], W5) and len(vac) == 2)
say(f'done ({time.time()-T0:.0f}s)')
if FAIL: say(f'CONTROLS FAILED {FAIL} -> exit 1'); sys.exit(1)
say('ALL ENGINE CONTROLS PASS'); sys.exit(0)
