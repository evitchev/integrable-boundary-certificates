"""q = 8 test, w3: the full kernel K_14 = ker ad I_3 at weight 14 mod d, on the complement basis C_14 (895 columns), by the Euler + fmpq_mat
route; GATES G2 (dim 2) and G3 (certified Sol-3 weight-14 eigenvalue at 9/7 in span eigen(K_14)); cross-checks X1-X5 (the recorded predictions).
Exit: 0 both gates pass; 1 a gate fails (INCONCLUSIVE) or I_3 is not one-dimensional."""
import os, sys, time, json
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import q8_lib as L
from q8_lib import K, A, B, say, PX, PI, PA, PB
import sympy as sp
OUT = sys.argv[1]   # output directory for w3_kernel.json (run.sh passes it)
say('w3_kernel.py sha256', L.self_hash(__file__))
T0 = time.time(); FAIL = []
def gate(name, ok):
    ok = bool(ok); say(f'GATE {name}: {ok}  ({time.time()-T0:.0f}s)')
    if not ok: FAIL.append(name)
def xcheck(name, ok): say(f'CROSS-CHECK {name}: {bool(ok)}  ({time.time()-T0:.0f}s)')
K.set_kappa(2)
c4 = K.reps_mod_d(K.commutant(4, [{A: F(1, 4), B: F(1, 4)}, {A: F(2), B: F(0)}, {A: F(0), B: F(2)}]), 4)
gate('I_3 dim 1 (recomputed)', len(c4) == 1); I3 = c4[0]
H7 = K.current(2, A); H7Z = K.current(2, B); G14 = K.mul(H7, H7Z)
C14 = L.complement_basis(14); say(f'   C_14: {len(C14)} monomials (p2(14) - p2(13) = {L.p2(14) - L.p2(13)})')
t1 = time.time(); vecs = []
for i, m in enumerate(C14):
    vecs.append(L.euler_vec(K.ope(I3, {m: F(1)}, poles={1}).get(1, {})))
    if (i + 1) % 100 == 0: say(f'   images + Euler: {i+1}/{len(C14)} ({time.time()-t1:.0f}s)')
keys = {k for v in vecs for k in v}
say(f'   images + Euler done: {time.time()-t1:.0f}s; matrix {len(keys)} x {len(C14)}, nonzeros {sum(len(v) for v in vecs)}')
t1 = time.time(); rel, rk = L.nullspace_Q(vecs); say(f'   exact nullspace over Q (fmpq_mat rref): rank {rk}, nullity {len(rel)} ({time.time()-t1:.0f}s)')
t1 = time.time(); rkp = [L.rank_mod(vecs, p) for p in (2147483647, 1000000007)]; say(f'   rank mod 2^31-1 and mod 10^9+7: {rkp} ({time.time()-t1:.0f}s)')
gate('G2 dim K_14 = 2 (VIR5d: kernel dimension 2 at t = 9/7)', len(rel) == 2)
xcheck('X5 nullity mod two primes equals the exact nullity', all(len(C14) - r == len(rel) for r in rkp))
reps = [{C14[j]: x[j] for j in range(len(C14)) if x[j] != 0} for x in rel]
say('   kernel representatives (C_14 coordinates): terms', [len(r) for r in reps])
t1 = time.time(); eK = [L.internal(L.eigen_fast(r), 2) for r in reps]; say(f'   eigenvalues of the kernel classes ({time.time()-t1:.0f}s):')
for i, e in enumerate(eK): say(f'   eigen(K_14 class {i}) = {e}')
cert = L.read_vev_profile(os.path.join(L.REPO, 'results/lab/anchor13/vev_profile_sol3_w14.json'), sp.Rational(9, 7))
say('   certified Sol-3 weight-14 eigenvalue at 9/7 (anchor13):', cert)
cs = sp.symbols('c0:2'); sol = sp.solve(sp.Poly(sp.expand(cert - sum(c * e for c, e in zip(cs, eK))), PX, PI).coeffs(), cs, dict=True) if len(eK) == 2 else []
gate('G3 the certified Sol-3 weight-14 eigenvalue at 9/7 lies in span eigen(K_14)', bool(sol))
if sol: say('   coefficients:', sol[0])
# X1: VIR5d's two classes at 9/7
kd = json.load(open(os.path.join(L.WS, 'inputs', 'z2_kernel_t9_7_K7.json')))
v5 = [{(int(k.split(',')[0]), int(k.split(',')[1])): sp.Rational(v) for k, v in c['eigenvalue'].items()} for c in kd['classes']]
mine = [sp.Poly(e, PX, PI).as_dict() for e in eK]
mons = sorted(set(m for d_ in v5 + mine for m in d_))
M = sp.Matrix([[d_.get(m, 0) for m in mons] for d_ in mine + v5])
xcheck(f'X1 span eigen(K_14) == span of VIR5d\'s two t = 9/7 classes (rank of the 4 x {len(mons)} matrix = {M.rank()}; mine {sp.Matrix(M[:2, :]).rank()}, VIR5d {sp.Matrix(M[2:, :]).rank()})', M.rank() == 2 and sp.Matrix(M[:2, :]).rank() == 2)
# X2: the antisymmetric class
eG = L.internal(L.eigen_fast(G14), 2)
if len(eK) == 2:
    E = cs[0] * eK[0] + cs[1] * eK[1]; anti = sp.solve(sp.Poly(sp.expand(L.swap_sym(E) + E), PX, PI).coeffs(), cs, dict=True)
    an = sp.expand(E.subs(anti[0])) if anti else 0
    free = [c for c in cs if c not in anti[0]] if anti else []
    say(f'   antisymmetric combinations: solution {anti}; free symbols {free}')
    xcheck('X2 exactly one antisymmetric class and it is proportional to eigen(G_14)', len(free) == 1 and L.proportional(an.subs(free[0], 1), eG, (PX, PI)))
s2 = sp.solve(sp.Poly(sp.expand(cert - cs[0] * eG), PX, PI).coeffs(), [cs[0]], dict=True)
xcheck('X3 the certified Sol-3 eigenvalue is NOT proportional to eigen(G_14) (expected: the Sol-3 class is symmetric)', not s2 and sp.expand(L.swap_sym(cert) - cert) == 0)
# X4: G_14's class in C_14 coordinates lies in the kernel span
t1 = time.time(); a = L.solve_unique_Q([L.euler_vec({c: F(1)}) for c in C14], L.euler_vec(G14))
if a is None: xcheck('X4 G_14 reduces to C_14 coordinates', False)
else:
    av = {j: a[j] for j in range(len(C14)) if a[j] != 0}; rv = [{j: x[j] for j in range(len(C14)) if x[j] != 0} for x in rel]
    xcheck(f'X4 G_14\'s class (C_14 coordinates, {len(av)} terms) lies in span K_14 (rank {L.rank_Q(rv + [av])} == {len(rv)}) ({time.time()-t1:.0f}s)', L.rank_Q(rv + [av]) == len(rv))
json.dump({'C14_size': len(C14), 'rank_Q': rk, 'rank_mod': rkp, 'dim_K14': len(rel), 'eigen_K14': [str(e) for e in eK], 'cert_w14_at_9_7': str(cert), 'eigen_G14': str(eG),
           'kernel_reps': [{str(k): str(v) for k, v in r.items()} for r in reps]}, open(os.path.join(OUT, 'w3_kernel.json'), 'w'), indent=1)
say(f'done ({time.time()-T0:.0f}s)')
if FAIL: say(f'GATE FAILED {FAIL} -> exit 1 (INCONCLUSIVE)'); sys.exit(1)
say('w3: gates G2, G3 pass'); sys.exit(0)
