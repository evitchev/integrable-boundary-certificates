"""level2.py -- B2 Level II: J-projections of the tower charges versus Solution 1's commutant in the coset frame.

Per fibre (t, s = p^2) and tower (Z or K):
  * charges from towers_<T>_s<..>_w*.json (w2 = 6, 8, 10, 12 as available), projected by coset.coset_project -> D(q) =
    sum_k q^k D_k, (u, v) densities with rational coefficients;
  * Solution 1's (u, v) commutant at the same weight (pure exponentials (1/2, -1/2), (1/2, 1/2), (-1, 0), metric (s, 4-s));
  * TEST (P-II): every D_k lies in that commutant mod d; at q = 0 the constant c_s with D_0 = c_s I_s^{cert} mod d
    (density level for I_3 via sol1_t<..>.json; eigenvalue level for I_5 via the certified kit eigenvalue);
  * commutation of D(q) with the pure exponential e^{-u} and with the projected W0 vertex, at several rational q;
  * info: v-parity (X-parity) of each D_0; S3, S4 (full frame) on the unprojected charge; (u,v)-components of S3, S4 on D_0;
    is Z's I_3 the same density as K's?
  * CONTROLS: wrong dictionary (tower at s' = s + 1/7 versus Solution 1 at s, compared in the (X, phi) frame for I_3 and
    at the eigenvalue level for I_5); perturbed screening vector (W0 momentum (-1, 0) -> (-1 - 1/7, 0)).
Usage: level2.py <tower> <t>
"""
import sys
sys.dont_write_bytecode = True
import glob
import json
import time
from fractions import Fraction as F

import sympy as sp

import gm
from gm import X, add, scale, clean
from super2_screening import Vertex, boson_vertex, basis, sector_candidates
import n2frame as N2
import coset
from towers import cur_from_json

tower, t = sys.argv[1], F(sys.argv[2])
s = 2 * (t - 1) / (t + 1)
tag = str(s).replace('/', '-')
ttag = str(t).replace('/', '-')
out = dict(tower=tower, t=str(t), s=str(s), cells={})
t0 = time.time()


def P(msg):
    print(msg, flush=True)


def load_charges(tw, sv):
    stag = str(sv).replace('/', '-')
    ch = {}
    for fn in glob.glob(f'towers_{tw}_s{stag}_w*.json'):
        if 'perturb' in fn:
            continue
        d = json.load(open(fn))
        for w2, r in d['results'].items():
            ch[int(w2)] = [cur_from_json(k) for k in r['kernel']]
    return ch


_GENS = {}


def sol1_gens(sv, w2, third=(F(-1), F(0))):
    """Solution 1's (u,v) commutant at doubled weight w2 (memoised); w2 <= 2: no genuine charges (weight-1 densities are
    total derivatives), returned as []."""
    key = (F(sv), w2, third)
    eta_uv = (F(sv), 4 - F(sv))
    if key not in _GENS:
        if w2 <= 2:
            _GENS[key] = []
        else:
            Vuv = [boson_vertex((F(1, 2), F(-1, 2))), boson_vertex((F(1, 2), F(1, 2))), boson_vertex(third)]
            cands = [{m: F(1)} for m in basis(w2, 2) if all(l[0] == 0 for l in m)]
            _GENS[key] = gm.joint_kernel(cands, Vuv, eta_uv)
    return _GENS[key], eta_uv


def split_weight(cur):
    """a projected q-coefficient is a sum of (u,v) densities of different doubled weights (the H~ contractions produce
    cylinder constants): split it into weight-homogeneous pieces {w2: current}."""
    out = {}
    for m, c in cur.items():
        w = gm.mono_weight2(m) if m else 0
        out.setdefault(w, {})[m] = c
    return out


def membership(Dq, w2, sv, third=(F(-1), F(0))):
    """every (q^k, weight-w) piece must lie in Solution 1's commutant at weight w (mod d); weight 0 (constants) and
    weight 2 (all total derivatives) are trivially fine.  Returns {k: {w: (member, dim)}}."""
    out = {}
    for k in sorted(Dq):
        out[k] = {}
        for w, cur in sorted(split_weight(Dq[k]).items()):
            if w <= 0:
                out[k][w] = (True, 'const')
            elif w == 2:
                out[k][w] = (gm.is_total_derivative(cur), 'weight-1')
            else:
                gens, _ = sol1_gens(sv, w, third)
                out[k][w] = (gm.in_span_mod_d(cur, gens), len(gens))
    return out


def all_members(mem):
    return all(v[0] for d in mem.values() for v in d.values())


def commute_orders(Dq, w2, Wq, b, eta_uv):
    """[oint sum_j q^j W_j e^{b.Y}, oint sum_k q^k D_k] = 0 order by order in q: at order n the residue
    sum_{j+k=n} Res(D_k, W_j e^{b.Y}) (one doubled weight) must be twisted-exact."""
    from super2_screening import is_twisted_exact, vertex_residue
    res = {}
    nmax = max(Dq) + max(Wq)
    for n in range(nmax + 1):
        byw = {}
        for j, W in Wq.items():
            k = n - j
            if k not in Dq or not Dq[k] or not W:
                continue
            V = Vertex(b, W)
            for w, piece in split_weight(Dq[k]).items():
                r = vertex_residue(piece, V, eta_uv)
                wt = w + V.w2W - 2
                byw[wt] = add((1, byw.get(wt, {})), (1, r))
        if not byw:
            continue
        ok = True
        for wt, R in byw.items():
            if wt <= 0:
                ok &= (not R)
            elif R:
                ok &= is_twisted_exact(R, b, wt, 2)
        res[n] = ok
    return res


def reduce_mod_d(cur, nf=2):
    r, piv, index = gm.mod_d_coords(cur, nf)
    return r


def proportionality(A, B, nf=2):
    """A = c B mod d ?  Exact linear solve over Q on the bosonic monomial basis: returns c (sympy Rational) or None."""
    if not A:
        return sp.Integer(0)
    w2 = gm.mono_weight2(next(iter(B)))
    mons = [m for m in basis(w2, nf) if all(l[0] == 0 for l in m)]
    lower = [m for m in basis(w2 - 2, nf) if all(l[0] == 0 for l in m)] if w2 > 2 else []
    idx = {m: i for i, m in enumerate(mons)}
    M = sp.zeros(len(mons), len(lower) + 1)
    for j, u in enumerate(lower):
        for mm, cc in gm.d_mono(u).items():
            M[idx[mm], j] = cc
    for mm, cc in B.items():
        M[idx[mm], len(lower)] = sp.Rational(cc)
    v = sp.zeros(len(mons), 1)
    for mm, cc in A.items():
        v[idx[mm], 0] = sp.Rational(cc)
    try:
        sol, params = M.gauss_jordan_solve(v)
    except ValueError:
        return None
    c = sol[len(lower), 0]
    if params:
        c = c.subs({pp: 0 for pp in params})
    return sp.nsimplify(c)


def vparity(cur):
    ps = {sum(1 for l in m if l[2] == 1) % 2 for m in cur}
    return ps


sv = F(s)
charges = load_charges(tower, sv)
sol1 = json.load(open(f'sol1_t{ttag}.json'))
I3_uv = cur_from_json(sol1['I3_uv'])
P(f"tower {tower}, t = {t}, s = {s}; charges at w2 = {sorted(charges)}; dims {[len(charges[w]) for w in sorted(charges)]}")
Bz = coset.Bosoniser(sv)
C = N2.currents(sv)
S = N2.kl_vertices(sv)
eta_full = N2.eta_of(sv)
PX, Ppi = sp.symbols('P_X pi')
p_ = sp.sqrt(sp.Rational(sv)); rho = p_ / 2 - 1 / p_; kappa = sp.sqrt(1 - sp.Rational(sv) / 4)
P2, Q2 = sp.symbols('P2 Q2')
dictP, dictQ = sol1['momentum_dictionary']
assert dictP == 'P2 = +P_X^2' and dictQ.startswith('Q2 = pi^2 - rho^2'), sol1['momentum_dictionary']
# evaluator momenta: P_X, and P_phi(evaluator) = pi; our P2 = P_X^2, Q2 = pi^2 - rho^2 (sol1.py, exact match of I_3)
subsPQ = {P2: PX ** 2, Q2: Ppi ** 2 - rho ** 2}
E_cert = {8: sp.expand(sp.sympify(sol1['E3_repo']).subs(subsPQ)), 12: sp.expand(sp.sympify(sol1['E5_alt']).subs(subsPQ))}
I5_top_alt = sp.Rational(sol1['I5_top_coeff_alt'])      # kit normalisation = alt / this (0 => kit undefined)
Pu, Pv = sp.symbols('P_u P_v')
# u = -p phi, v = -2 kappa X on the fields, hence on the evaluator momenta (same recipe, linear map): P_u = -p pi, P_v = -2 kappa P_X
mom_map = {Pu: -p_ * Ppi, Pv: -2 * kappa * PX}


def eig_uv(D, eta_uv):
    return sp.expand(gm.vev_ns(D, eta_uv, [Pu, Pv]).subs(mom_map))


def eig_ratio(E, Ecert):
    """E = c Ecert as polynomials in (P_X, pi) with sqrt coefficients: returns c or None."""
    c = sp.Symbol('c')
    eqs = sp.Poly(sp.expand(E - c * Ecert), PX, Ppi).coeffs()
    sol = sp.solve(eqs, c, dict=True)
    return sol[0][c] if sol else None


for w2 in sorted(charges):
    if w2 < 6:
        continue
    for i, D in enumerate(charges[w2]):
        cell = {}
        Dq, info = coset.coset_project(D, sv, Bz)
        gens, eta_uv = sol1_gens(sv, w2)
        qmax = max(Dq)
        member = membership(Dq, w2, sv)
        P(f"\n== w2 = {w2} (spin {w2 // 2 - 1}) charge {i}: J-charge sectors (size, total derivative?) {info}; "
          f"projected q-powers {sorted(Dq)}; Solution 1 (u,v) commutant dim {len(gens)}")
        P(f"   (q^k, weight-w) pieces in Solution 1's commutant at weight w mod d [k: {{w: (member, dim)}}]: {member}; ALL: {all_members(member)}")
        member = {k: all(v[0] for v in d.values()) for k, d in member.items()}
        cell.update(dim_sol1=len(gens), q_powers=sorted(Dq), membership={str(k): v for k, v in member.items()},
                    nonneutral={str(k): v for k, v in info.items()})
        D0full = Dq.get(0, {})
        D0 = split_weight(D0full).get(w2, {}) if D0full else {}     # the top (weight-w2) piece of the q = 0 part
        cell['D0_weights'] = sorted(split_weight(D0full)) if D0full else []
        P(f"   q = 0 part has weight pieces {cell['D0_weights']}; top piece (weight {w2}) has {len(D0)} monomials")
        cell['vparity_D0'] = sorted(vparity(D0)) if D0 else None
        cell['vparity_gens'] = [sorted(vparity(g)) for g in gens]
        P(f"   v-parity (= X-parity) of D_0: {cell['vparity_D0']}; of the commutant generators: {cell['vparity_gens']}")
        if not D0:
            P("   D_0 = 0 (the q = 0 projection vanishes)")
        # constants
        if w2 == 8 and D0:
            c3 = proportionality(D0, I3_uv)
            cell['c3_density'] = str(c3)
            P(f"   q = 0: D_0 = c_3 * I_3^cert mod d with c_3 = {c3}")
            # the q-dependence: D_k = c_k I_3 ?
            cs = {}
            for k in sorted(Dq):
                for w, piece in sorted(split_weight(Dq[k]).items()):
                    if w == 8:
                        cs[(k, w)] = ('I_3^cert', proportionality(piece, I3_uv))
                    elif w > 2:
                        g, _ = sol1_gens(sv, w)
                        cs[(k, w)] = (f'gen w2={w}', proportionality(piece, g[0]) if len(g) == 1 else 'dim != 1')
                    else:
                        cs[(k, w)] = ('trivial', str(piece))
            cell['ck'] = {str(k): str(v) for k, v in cs.items()}
            P(f"   (q^k, weight w) piece = c x (generator at that weight): {cs}")
        if w2 in (8, 12) and D0:
            E = eig_uv(D0, eta_uv)
            ratio = eig_ratio(E, E_cert[w2])
            cell['c_eigen'] = str(ratio)
            nrm = 'kit (dX)^4 normalisation of cyl_P4' if w2 == 8 else 'ALT normalisation (P2^2 Q2 coefficient = 1)'
            P(f"   eigenvalue level: vev(D_0)(P_u = -p pi, P_v = -2 kappa P_X) = c * E_cert [{nrm}]: c = {ratio}")
            if w2 == 12 and ratio is not None:
                ckit = (ratio * I5_top_alt) if I5_top_alt != 0 else None
                cell['c5_kit'] = str(ckit)
                P(f"   c_5 in the kit normalisation ((dX)^6 = 1): {ckit if ckit is not None else 'undefined at this fibre ((dX)^6 coefficient = 0)'}")
            if len(gens) == 1 and ratio is not None:
                cell['c_from_eigen_and_1dim'] = str(ratio)
        if w2 == 6 and D0 and gens:
            cs = {}
            for k in sorted(Dq):
                for w, piece in sorted(split_weight(Dq[k]).items()):
                    g, _ = sol1_gens(sv, w) if w > 2 else ([], None)
                    cs[(k, w)] = proportionality(piece, g[0]) if (w > 2 and len(g) == 1) else ('trivial' if w <= 2 else 'dim != 1')
            P(f"   spin 2: (q^k, weight w) piece = c x (Solution 1 generator at weight w) mod d: {cs}")
            cell['ck_vs_gen'] = {str(k): str(v) for k, v in cs.items()}
        # commutation with e^{-u} and with the projected W0 at rational q
        W0p = coset.project_vertex_prefactor(N2.w0_vertex(sv, 'G')[1].W, sv, Bz)
        one = {(): F(1)}
        comm = dict(pure=commute_orders(Dq, w2, {0: one}, (F(-1), F(0)), eta_uv),
                    W0proj=commute_orders(Dq, w2, W0p, (F(-1), F(0)), eta_uv),
                    S1=commute_orders(Dq, w2, {0: one}, (F(1, 2), F(-1, 2)), eta_uv),
                    S2=commute_orders(Dq, w2, {0: one}, (F(1, 2), F(1, 2)), eta_uv))
        cell['commutation_by_order'] = {k: {str(n): v for n, v in d.items()} for k, d in comm.items()}
        P(f"   [oint V, oint D(q)] = 0 order by order in q for V = e^(-u) pure, W0 projected (prefactor {W0p}), S1uv, S2uv: {comm}")
        # info: S3, S4 on the unprojected charge, and their (u,v)-components on D_0
        info34 = {k: gm.screens(D, S[k], eta_full) for k in (3, 4)}
        s34uv = {3: boson_vertex((F(-1, 2), (sv + 4) / (2 * (sv - 4)))), 4: boson_vertex((F(-1, 2), -(sv + 4) / (2 * (sv - 4))))}
        info34uv = {k: gm.screens(D0, v, eta_uv) if D0 else None for k, v in s34uv.items()}
        cell['S3S4_full'] = info34
        cell['S3S4_uv_on_D0'] = info34uv
        P(f"   info: unprojected charge commutes with S3, S4: {info34}; (u,v)-components of S3, S4 on D_0: {info34uv}")
        out['cells'][f'w{w2}_{i}'] = cell

# Z vs K at spin 3 and 5: same density?
if tower == 'Z':
    chK = load_charges('K', sv)
    for w2 in (8, 12):
        if w2 in chK and w2 in charges and chK[w2] and charges[w2]:
            same = gm.in_span_mod_d(charges[w2][0], chK[w2])
            P(f"\ninfo: Z's charge at w2 = {w2} equals K's mod d: {same}")
            out[f'Z_equals_K_w{w2}'] = same

# ---------------------------------------------------------------- controls
P("\n== CONTROLS")
# (i) wrong dictionary: tower at s' = s + 1/7
sp_ = sv + F(1, 7)
Vq = N2.kl_vertices(sp_)
Vl = [Vq[1], Vq[2], Vq[3], Vq[4]] if tower == 'K' else [Vq[1], Vq[2], N2.w0_vertex(sp_, 'G')[1]]
Bz2 = coset.Bosoniser(sp_)
ctrl = {}
I3_Xphi = {tuple(tuple(l) for l in m): sp.sympify(c) for m, c in sol1['I3_Xphi']}
for w2 in (8, 12):
    ker = gm.joint_kernel(sector_candidates(w2, 'all'), Vl, N2.eta_of(sp_))
    if not ker:
        ctrl[f'w{w2}'] = 'no charge'
        continue
    Dq2, _ = coset.coset_project(ker[0], sp_, Bz2)
    D02 = split_weight(Dq2.get(0, {})).get(w2, {})
    gens_s, eta_uv_s = sol1_gens(sv, w2)
    # (a) in the (u,v) frame of s' the density is rational; Solution 1's commutant at s is in the (u,v) frame of s:
    #     compare in the (X, phi) frame: u' = -p' phi, v' = -2 kappa' X
    pp = sp.sqrt(sp.Rational(sp_)); kp = sp.sqrt(1 - sp.Rational(sp_) / 4)
    DX = {}
    for m, c in D02.items():
        nu = sum(1 for l in m if l[2] == 0); nv = sum(1 for l in m if l[2] == 1)
        key = tuple(sorted((0, l[1], 1 if l[2] == 0 else 0) for l in m))   # u -> phi (index 1), v -> X (index 0)
        DX[key] = DX.get(key, 0) + sp.Rational(c) * (-pp) ** nu * (-2 * kp) ** nv
    DX = {m: sp.expand(c) for m, c in DX.items() if sp.expand(c) != 0}
    if w2 == 8:
        # membership test of DX in span(I3_Xphi) + Im d over the algebraic numbers
        mons = sorted(set(DX) | set(I3_Xphi))
        w2x = gm.mono_weight2(next(iter(I3_Xphi)))
        lower = [m for m in basis(w2x - 2, 2) if all(l[0] == 0 for l in m)]
        cols = []
        for u in lower:
            dm = {}
            for mm, cc in gm.d_mono(u).items():
                dm[mm] = dm.get(mm, 0) + cc
            cols.append(dm)
        allm = sorted(set(mons) | {mm for dm in cols for mm in dm})
        idx = {m: i for i, m in enumerate(allm)}
        M = sp.zeros(len(allm), len(cols) + 1)
        for j, dm in enumerate(cols):
            for mm, cc in dm.items():
                M[idx[mm], j] = cc
        for mm, cc in I3_Xphi.items():
            M[idx[mm], len(cols)] = cc
        v = sp.zeros(len(allm), 1)
        for mm, cc in DX.items():
            v[idx[mm], 0] = cc
        rk0 = M.rank(simplify=True)
        rk1 = M.row_join(v).rank(simplify=True)
        inspan = (rk0 == rk1)
        # the same test for the REAL projection at s (must be True)
        ctrl[f'w{w2}_wrong_dictionary_in_span'] = inspan
        P(f"   control wrong dictionary (s' = {sp_}) at w2 = 8: D_0(s') in span(I_3^cert(t)) + Im d in the (X, phi) frame: {inspan}  (must be False)")
    else:
        E2 = sp.expand(gm.vev_ns(D02, (sp_, 4 - sp_), [Pu, Pv]).subs({Pu: -pp * Ppi, Pv: -2 * kp * PX}))
        r = eig_ratio(E2, E_cert[12])
        ctrl['w12_wrong_dictionary_eigen_ratio'] = str(r)
        P(f"   control wrong dictionary at w2 = 12: eigenvalue proportional to the certified I_5: {r}  (must be None)")
# also the real projection at s in the (X, phi) frame for w2 = 8 (consistency with the rational test)
Dq, _ = coset.coset_project(charges[8][0], sv, Bz)
D0 = split_weight(Dq[0]).get(8, {})
DXr = {}
for m, c in D0.items():
    nu = sum(1 for l in m if l[2] == 0); nv = sum(1 for l in m if l[2] == 1)
    key = tuple(sorted((0, l[1], 1 if l[2] == 0 else 0) for l in m))
    DXr[key] = DXr.get(key, 0) + sp.Rational(c) * (-p_) ** nu * (-2 * kappa) ** nv
DXr = {m: sp.expand(c) for m, c in DXr.items() if sp.expand(c) != 0}
mons = sorted(set(DXr) | set(I3_Xphi))
lower = [m for m in basis(6, 2) if all(l[0] == 0 for l in m)]
cols = [dict(gm.d_mono(u)) for u in lower]
allm = sorted(set(mons) | {mm for dm in cols for mm in dm})
idx = {m: i for i, m in enumerate(allm)}
M = sp.zeros(len(allm), len(cols) + 1)
for j, dm in enumerate(cols):
    for mm, cc in dm.items():
        M[idx[mm], j] = cc
for mm, cc in I3_Xphi.items():
    M[idx[mm], len(cols)] = cc
v = sp.zeros(len(allm), 1)
for mm, cc in DXr.items():
    v[idx[mm], 0] = cc
inspan_real = (M.rank(simplify=True) == M.row_join(v).rank(simplify=True))
P(f"   real projection at s in the (X, phi) frame: D_0 in span(I_3^cert) + Im d: {inspan_real}  (must be True)")
ctrl['w8_real_in_span_Xphi'] = inspan_real
# (ii) perturbed screening vector (last vertex's momentum shifted by -1/7 along u)
Vp = [Vq_ for Vq_ in (N2.kl_vertices(sv)[1], N2.kl_vertices(sv)[2])]
if tower == 'Z':
    w0 = N2.w0_vertex(sv, 'G')[1]
    Vp.append(Vertex((F(-1) - F(1, 7), F(0)), w0.W, label='W0 perturbed'))
else:
    S_ = N2.kl_vertices(sv)
    Vp += [S_[3], Vertex((S_[4].b[0] - F(1, 7), S_[4].b[1]), S_[4].W, label='S4 perturbed')]
for w2 in (8, 12):
    ker = gm.joint_kernel(sector_candidates(w2, 'all'), Vp, eta_full)
    if not ker:
        ctrl[f'w{w2}_perturbed_screening'] = 'no charge (commutant empty) -> fails'
        P(f"   control perturbed screening at w2 = {w2}: commutant dim 0 (no charge to project) -> Level II fails")
        continue
    Dqp, _ = coset.coset_project(ker[0], sv, Bz)
    memf = membership(Dqp, w2, sv)
    mem = {k: all(v[0] for v in d.values()) for k, d in memf.items()}
    ctrl[f'w{w2}_perturbed_screening'] = {str(k): v for k, v in mem.items()}
    P(f"   control perturbed screening at w2 = {w2}: commutant dim {len(ker)}; membership {mem}  (must contain False)")
# (iii) perturbed Solution 1 vector in the (u,v) comparison
for w2 in (8, 12):
    third = (F(-1) - F(1, 7), F(0))
    gensp, _ = sol1_gens(sv, w2, third)
    Dq, _ = coset.coset_project(charges[w2][0], sv, Bz)
    memf = membership(Dq, w2, sv, third)
    mem = {k: all(v[0] for v in d.values()) for k, d in memf.items()}
    ctrl[f'w{w2}_perturbed_sol1_vector'] = str(mem)
    P(f"   control perturbed Solution 1 vector at w2 = {w2}: commutant dim {len(gensp)}; membership {mem}  (must fail)")
out['controls'] = ctrl
json.dump(out, open(f'level2_{tower}_t{ttag}.json', 'w'), indent=1, default=str)
P(f"done in {time.time() - t0:.1f}s")
