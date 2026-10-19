"""sol1.py -- B2: the certified Solution 1 data in the coset frame.

(1) I_3: code/mixed_engine.cyl_P4(1, N, sc) at the fibre t (N = (t^2-25)/(t^2-1), sc = -24t/(t^2-1)), a density in
    (X-derivative orders, Y-pairs).  Y-pairs -> Vir_N: (11) = 2 T_Y, and (11)(11) - 2(22) = 4 (T_Y T_Y) mod d (point-split
    composite; re-verified here with the VIR1 Y-pair OPE engine results/lab/large_n/vir1/vir_lib.py, identical copy).  Then T_Y = 1/2 (d phi)^2 +
    rho d^2 phi, rho = p/2 - 1/p, metric (1, 1) for (X, phi), composites by gm.normal_prod.  Then the isometry of Level I:
    u = -p phi, v = -2 kappa X (kappa^2 = 1 - p^2/4): d^m phi = -d^m u / p, d^m X = -d^m v/(2 kappa); the result must be
    RATIONAL (asserted).
(2) The (u, v) commutant of Solution 1's screening vectors (1/2, -1/2), (1/2, 1/2), (-1, 0) [pure exponentials, metric
    (s, 4-s)] at doubled weights 6, 8, 12 (bosonic candidates): dimensions and generators; I_3^{uv} membership.
(3) Eigenvalues: the certified vacuum eigenvalue of I_3 (repo cyl_vev.vev_density on cyl_P4) compared with the (X, phi)
    density's NS cylinder eigenvalue (gm.vev_ns) to pin the momentum dictionary P2, Q2 <-> (P_X, P_phi); the certified I_5
    (kernel of ad I_3 at weight 6, normalised (dX)^6 = 1 as the kit) and its eigenvalue versus the kit table
    results/lab/vev/vev_sol1_w6.json at the fibre.
Usage: sol1.py <t>   (e.g. 2, 7/2, 11/3, 13/5)
"""
import sys
sys.dont_write_bytecode = True
import os as _os_kit; _KITROOT = _os_kit.path.join(_os_kit.path.dirname(_os_kit.path.abspath(__file__)), '..', '..', '..', '..', '..') + '/'; sys.path.insert(0, _KITROOT + 'results/lab/large_n/vir1')  # KIT PATCH: was an absolute path to a private copy of vir1/ (vir_lib.py identical, sha256 fc851415...)
import json
from fractions import Fraction as F

import sympy as sp

import gm
from gm import X, add, scale, clean
from super2_screening import boson_vertex, basis
import mixed_engine as ME
import cyl_vev as CV
import vir_lib as VL

t = F(sys.argv[1])
N = (t * t - 25) / (t * t - 1)
sc = F(-24) * t / (t * t - 1)
s = 2 * (t - 1) / (t + 1)            # p^2
p = sp.sqrt(sp.Rational(s))
rho = p / 2 - 1 / p
kappa = sp.sqrt(1 - sp.Rational(s) / 4)
print(f"t = {t}: N = {N}, sc = {sc}, p^2 = s = {s}, rho = {rho}, kappa^2 = {1 - s/4}, c_{{X,phi}} = {2 - 12*sp.expand(rho**2)}")
out = dict(t=str(t), N=str(N), s=str(s))

# ---------------------------------------------------------------- (1) certified I_3 -> (X, phi) -> (u, v)
I3m = ME.cyl_P4(1, N, sc)
print("certified I_3 (mixed basis):", I3m)
expect_keys = {ME._M((), [(1, 1), (1, 1)]), ME._M((), [(2, 2)]), ME._M((2, 2)), ME._M((1, 1), [(1, 1)]), ME._M((1, 1, 1, 1))}
assert set(I3m) == expect_keys, set(I3m)
assert I3m[ME._M((), [(1, 1), (1, 1)])] == 1 and I3m[ME._M((), [(2, 2)])] == -2
c4, c22, cxy = I3m[ME._M((1, 1, 1, 1))], I3m[ME._M((2, 2))], I3m[ME._M((1, 1), [(1, 1)])]
# re-verify (11)(11) - 2(22) == 4 (T_Y T_Y) mod d in the Y-pair algebra (VIR1 engine, loop variable = N)
TY = {VL.canon((), [(1, 1)]): F(1, 2)}
TT = VL.ope_vec(TY, TY, 0, N)                      # regular term = point-split (T_Y T_Y)
lhs = {VL.canon((), [(1, 1), (1, 1)]): F(1), VL.canon((), [(2, 2)]): F(-2)}
diff = dict(lhs)
for m, c in TT.items():
    diff[m] = diff.get(m, 0) - 4 * c
diff = {m: c for m, c in diff.items() if c}
# membership in Im d at weight 4 (Y-only): d of weight-3 Y-only monomials
B3 = [m for m in VL.basis(3) if not m[0]]
cols = [VL.d_vec({m: F(1)}) for m in B3]
B4 = VL.basis(4)
idx = {m: i for i, m in enumerate(B4)}
rows = [[F(0)] * len(B4) for _ in cols]
for r, col in enumerate(cols):
    for m, c in col.items():
        rows[r][idx[m]] = c
target = [F(0)] * len(B4)
for m, c in diff.items():
    target[idx[m]] = c
rk0 = VL.rank(rows)
rk1 = VL.rank(rows + [target])
print(f"(11)(11) - 2(22) - 4(T_Y T_Y) is a total derivative in the Y-pair algebra: {rk0 == rk1}")
assert rk0 == rk1
# (X, phi) realisation: field 0 = X, field 1 = phi, metric (1, 1)
eta2 = (1, 1)
TYphi = {(X(1, 1), X(1, 1)): sp.Rational(1, 2), (X(2, 1),): rho}
TTphi = gm.normal_prod(TYphi, TYphi, eta2)
dX2 = {(X(1, 0), X(1, 0)): sp.Integer(1)}
I3_Xphi = add((4, TTphi), (cxy * 2, gm.normal_prod(dX2, TYphi, eta2)),
              (sp.Rational(c4), {(X(1, 0),) * 4: sp.Integer(1)}), (sp.Rational(c22), {(X(2, 0), X(2, 0)): sp.Integer(1)}))
I3_Xphi = {m: sp.expand(c) for m, c in I3_Xphi.items()}
print("I_3 in (X, phi) letters:", {m: str(c) for m, c in I3_Xphi.items()})
# (u, v): X = -v/(2 kappa), phi = -u/p
I3_uv = {}
for m, c in I3_Xphi.items():
    nX = sum(1 for l in m if l[2] == 0)
    nphi = sum(1 for l in m if l[2] == 1)
    cc = sp.expand(c * (-1 / (2 * kappa)) ** nX * (-1 / p) ** nphi)
    cc = sp.nsimplify(sp.simplify(cc))
    assert cc.is_Rational, (m, cc)
    key = tuple(sorted((0, l[1], 1 if l[2] == 0 else 0) for l in m))     # X -> v (index 1), phi -> u (index 0)
    I3_uv[key] = I3_uv.get(key, 0) + F(int(cc.p), int(cc.q))
I3_uv = clean(I3_uv)
print("I_3 in (u, v) letters (rational):", I3_uv)
out['I3_uv'] = [[list(map(list, m)), str(c)] for m, c in sorted(I3_uv.items())]
out['I3_Xphi'] = [[list(map(list, m)), str(c)] for m, c in sorted(I3_Xphi.items())]

# ---------------------------------------------------------------- (2) (u, v) commutant of Solution 1's vectors
eta_uv = (F(s), 4 - F(s))
Vuv = [boson_vertex((F(1, 2), F(-1, 2)), label='S1uv'), boson_vertex((F(1, 2), F(1, 2)), label='S2uv'),
       boson_vertex((F(-1), F(0)), label='W0uv')]
gens = {}
for w2 in (4, 6, 8, 10, 12):
    cands = [{m: F(1)} for m in basis(w2, 2) if all(l[0] == 0 for l in m)]
    ker = gm.joint_kernel(cands, Vuv, eta_uv)
    gens[w2] = ker
    xodd = [sum(1 for l in next(iter(k)) if l[2] == 1) % 2 for k in ker]
    print(f"(u,v) commutant of Solution 1's vectors at w2 = {w2}: dim {len(ker)}; v-parity of generators {xodd}")
out['uv_commutant_dims'] = {w2: len(k) for w2, k in gens.items()}
out['uv_commutant'] = {w2: [[[list(map(list, m)), str(c)] for m, c in sorted(k.items())] for k in ker] for w2, ker in gens.items()}
print("certified I_3 (u,v) in the w2 = 8 commutant mod d:", gm.in_span_mod_d(I3_uv, gens[8]))
assert gm.in_span_mod_d(I3_uv, gens[8])

# ---------------------------------------------------------------- (3) eigenvalues and the momentum dictionary
P2, Q2 = sp.symbols('P2 Q2')
E3_repo = sp.expand(CV.vev_density(I3m, sp.Rational(N), P2, Q2))
PX, Pphi = sp.symbols('P_X P_phi')
E3_mine = gm.vev_ns(I3_Xphi, eta2, [PX, Pphi])
print("repo eigenvalue of certified I_3 (P2, Q2):", E3_repo)
print("my (X, phi) eigenvalue (P_X, P_phi):", sp.expand(E3_mine))
dict_found = None
for nameP, sP in (('P2 = +P_X^2', 1), ('P2 = -P_X^2', -1)):
    for nameQ, Qexpr in (('Q2 = pi^2 - rho^2 with pi = P_phi(evaluator)', Pphi ** 2 - rho ** 2),
                         ('Q2 = rho^2 - pi^2 with pi = P_phi(evaluator)', rho ** 2 - Pphi ** 2),
                         ('Q2 = (P_phi - rho)^2 - rho^2', (Pphi - rho) ** 2 - rho ** 2), ('Q2 = rho^2 - (P_phi - rho)^2', rho ** 2 - (Pphi - rho) ** 2)):
        cand = sp.expand(E3_repo.subs({P2: sP * PX ** 2, Q2: Qexpr}))
        if sp.simplify(cand - E3_mine) == 0:
            dict_found = (nameP, nameQ)
            print(f"MOMENTUM DICTIONARY: {nameP}, {nameQ}, pi = P_phi - rho  (exact match)")
print("dictionary:", dict_found)
out['momentum_dictionary'] = dict_found
out['E3_repo'] = str(E3_repo)
# certified I_5: kernel of ad I_3 at weight 6, kit normalisation (dX)^6 = 1
K6 = ME.genuine_kernel_mixed(I3m, ME.gen_mixed_basis(6), N)
print(f"ad I_3 kernel at weight 6: dim {len(K6)}")
assert len(K6) == 1
I5m = K6[0]
top = ME._M((1, 1, 1, 1, 1, 1))
# ALT normalisation (defined at every fibre): vev coefficient of Q2^3 = 1.  Kit normalisation ((dX)^6 = 1) where the
# (dX)^6 coefficient is non-zero (it vanishes at t = 11/3: the kit's rational functions carry (3t - 11) denominators).
E5_raw = sp.expand(CV.vev_density(I5m, sp.Rational(N), P2, Q2))
k_alt = sp.Poly(E5_raw, P2, Q2).coeff_monomial(Q2 ** 3)      # non-zero at every fibre of the record (kit: -3(t-1)^2/((t-17)(3t-11)))
assert k_alt != 0
I5m = {m: F(c) / F(int(k_alt.p), int(k_alt.q)) for m, c in I5m.items()}
E5_alt = sp.expand(E5_raw / k_alt)
top_coeff = I5m.get(top, F(0))
print(f"certified I_5: (dX)^6 coefficient in the ALT normalisation (Q2^3 coefficient = 1): {top_coeff}; kit normalisation defined: {top_coeff != 0}")
out['I5_top_coeff_alt'] = str(top_coeff)
if top_coeff != 0:
    I5kit = {m: c / top_coeff for m, c in I5m.items()}
    E5_repo = sp.expand(E5_alt / sp.Rational(top_coeff))
else:
    I5kit = None
    E5_repo = None
out['E5_alt'] = str(E5_alt)
if E5_repo is None:
    print('kit normalisation singular at this fibre; E5_alt recorded; kit comparison skipped')
    out['E5_repo'] = None
    out['E5_kit_match'] = None
    out['I5_mixed_alt'] = [[str(m), str(c)] for m, c in I5m.items()]
    json.dump(out, open(f"sol1_t{str(t).replace('/', '-')}.json", 'w'), indent=1)
    print('written')
    sys.exit(0)
kit = json.load(open(_KITROOT + 'results/lab/vev/vev_sol1_w6.json'))  # KIT PATCH: repo root from __file__ (was an absolute path)
tt = sp.Symbol('t')
E5_kit = sum(sp.sympify(v).subs(tt, sp.Rational(t)) * sp.Symbol('P') ** int(k.split()[0][2:]) * sp.Symbol('Q') ** int(k.split()[1][2:])
             for k, v in kit['coefficients'].items())
# kit: article momenta P_art = iP -> P_art^2 = -P2 (our), Q_art^2 = -Q2 ; compare
Pa, Qa = sp.Symbol('P'), sp.Symbol('Q')
E5_kit_ours = sp.expand(E5_kit.subs({Pa: sp.I * sp.sqrt(P2), Qa: sp.I * sp.sqrt(Q2)}))
match = sp.simplify(E5_kit_ours - E5_repo) == 0
print(f"certified I_5 eigenvalue (repo recipe) == kit table at t = {t} (article P = i P_ours): {match}")
if not match:
    E5_kit_ours2 = sp.expand(E5_kit.subs({Pa: sp.sqrt(P2), Qa: sp.sqrt(Q2)}))
    print("   with P_art = P_ours instead:", sp.simplify(E5_kit_ours2 - E5_repo) == 0)
out['E3_repo'] = str(E3_repo)
out['E5_repo'] = str(E5_repo)
out['E5_kit_match'] = bool(match)
out['I5_mixed_alt'] = [[str(m), str(c)] for m, c in I5m.items()]
json.dump(out, open(f"sol1_t{str(t).replace('/', '-')}.json", 'w'), indent=1)
print('written')
