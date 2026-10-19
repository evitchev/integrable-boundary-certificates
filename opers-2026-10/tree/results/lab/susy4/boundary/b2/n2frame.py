"""n2frame.py -- B2: the N = 2 free-field realisation in the seat's frame, V1, vertices.

FRAME.  Two N = 1 superfields (X_j, psi_j), j = 0, 1, metric eta = (s, -s), s = p^2 (symbolic or rational).
Null screening momenta m_1 = (1/2, 1/2), m_2 = (1/2, -1/2) (m_i.m_i = 0, m_1.m_2 = s/2), m_3 = -m_1, m_4 = -m_2.
T = (1/s)[1/2 (dX_0)^2 - 1/2 (dX_1)^2 - 1/2 psi_0 dpsi_0 + 1/2 psi_1 dpsi_1] + lam_0 d^2 X_0 with lam_0 fixed by
h(SUSY vertex at m_1) = h(at m_2) = 1  (computed, not assumed).  G = (1/s)[psi_0 dX_0 - psi_1 dX_1] + 2 lam_0 dpsi_0.
G' (second supercurrent) and J are DERIVED by an ansatz solved against the N = 2 OPEs (symbolic s).

Tasks: v1 (symbolic-s OPE closure, c), weights (symbolic), vertices (commutation at rational s), w0 (Zeitlin's third
vertex as Res(G, bottom), its weight and its commutation with G_0, G'_0, J_0).
"""
import sys
sys.dont_write_bytecode = True
import json
import time
from fractions import Fraction as F

import sympy as sp

import gm
from gm import X, PSI, add, scale, d_cur, clean, residue, normal_prod, ope_all, screens
from super2_screening import Vertex, vertex_residue, susy_vertex, fermion_vertex, boson_vertex

s = sp.Symbol('s', positive=True)   # s = p^2


def eta_of(sv):
    return (sv, -sv)


def currents(sv, lam0=None, derive=False):
    """T, G, G', J at metric (sv, -sv).  lam0 default -1/sv (the weight-1 condition, verified in task weights)."""
    eta = eta_of(sv)
    if lam0 is None:
        lam0 = -1 / sv if not isinstance(sv, F) else F(-1) / sv
    lam = [lam0, 0]
    T = gm.T_eta(eta, lam)
    G = gm.G_eta(eta, lam)
    # g' := G'/i (real in this indefinite frame; derived by derive_gprime_j at symbolic s, asserted in task v1):
    #     g' = (1/s)[-psi_0 dX_1 + psi_1 dX_0 - 2 dpsi_1]
    # complex supercurrents G+- = (G -+ g')/sqrt2 are REAL; Gp2 = G - g', Gm2 = G + g' (= sqrt2 G+-).
    # J := (Gp2 Gm2)_2 / 4  (so that G+ G- ~ (2c/3)/z^3 + 2J/z^2 + (2T + dJ)/z), computed on the fly.
    inv = 1 / sv if not isinstance(sv, F) else F(1) / sv
    gp = clean({(X(1, 1), PSI(0, 0)): -inv, (X(1, 0), PSI(0, 1)): inv, (PSI(1, 1),): -2 * inv})
    Gp2 = add((1, G), (-1, gp))
    Gm2 = add((1, G), (1, gp))
    J = scale(gm.ope_coeff(Gp2, Gm2, 2, eta), F(1, 4) if isinstance(sv, F) else sp.Rational(1, 4))
    return dict(T=T, G=G, gp=gp, Gp2=Gp2, Gm2=Gm2, J=J, eta=eta, lam0=lam0)


def same(a, b):
    return not add((1, a), (-1, b))


def derive_gprime_j():
    """Ansatz G' = (1/s)[a1 psi_0 dX_1 + a2 psi_1 dX_0 + a3 dpsi_1 + a4 dpsi_0]; impose G'G' = GG (all poles) and
    GG' with no cubic pole and simple pole = 1/2 d(double pole).  Returns the solutions and J := (double pole)/2."""
    a1, a2, a3, a4 = sp.symbols('a1 a2 a3 a4')
    eta = eta_of(s)
    C = currents(s)
    G, T = C['G'], C['T']
    Gp = {(X(1, 1), PSI(0, 0)): a1 / s, (X(1, 0), PSI(0, 1)): a2 / s, (PSI(1, 1),): a3 / s, (PSI(1, 0),): a4 / s}
    GG = ope_all(G, G, eta, 4)
    GpGp = ope_all(Gp, Gp, eta, 4)
    GGp = ope_all(G, Gp, eta, 4)
    eqs = []
    for n in range(1, 5):
        for m in set(GG[n]) | set(GpGp[n]):
            eqs.append(sp.expand(GpGp[n].get(m, 0) - GG[n].get(m, 0)))
    for m in GGp[3]:
        eqs.append(sp.expand(GGp[3][m]))
    for m in GGp[4]:
        eqs.append(sp.expand(GGp[4][m]))
    dJ2 = d_cur(GGp[2])
    for m in set(GGp[1]) | set(dJ2):
        eqs.append(sp.expand(GGp[1].get(m, 0) - sp.Rational(1, 2) * dJ2.get(m, 0)))
    eqs = [e for e in eqs if e != 0]
    sols = sp.solve(eqs, [a1, a2, a3, a4], dict=True)
    out = []
    for sol in sols:
        gps = clean({m: sp.simplify(c.subs(sol) / sp.I) for m, c in Gp.items()})
        out.append((sol, gps))
    return out


def task_v1(log):
    t0 = time.time()
    eta = eta_of(s)
    sols = derive_gprime_j()
    log(f"derive G': {len(sols)} solutions of the ansatz (G' = i g'): " + '; '.join(str(x[0]) for x in sols))
    C = currents(s)
    T, G, gp, Gp2, Gm2, J = C['T'], C['G'], C['gp'], C['Gp2'], C['Gm2'], C['J']
    found = [k for k, (sol, gps) in enumerate(sols) if same(gps, gp)]
    log(f"the fixed real g' of currents() is solution index {found} (the other is g' -> -g')")
    assert found, "fixed g' not among the derived solutions"
    log(f"  g' = {gp}")
    log(f"  J  = (G+ G-)_2 / 2 in sqrt2-normalisation = {J}")
    c = sp.expand(3 - 12 / s)
    one = {(): sp.Integer(1)}
    res = {}
    allok = True

    def chk(name, got, want):
        nonlocal allok
        ok = all(sp.simplify(got.get(m, 0) - want.get(m, 0)) == 0 for m in set(got) | set(want))
        res[name] = ok
        allok &= ok
        log(f"  V1 {name}: {'OK' if ok else 'FAIL'}" + ('' if ok else f'  got {got} want {want}'))

    def poles(A, B):
        return ope_all(A, B, eta, 5)

    TT = poles(T, T)
    chk('TT/4 = c/2', TT[4], scale(one, c / 2)); chk('TT/3 = 0', TT[3], {}); chk('TT/2 = 2T', TT[2], scale(T, 2))
    chk('TT/1 = dT', TT[1], d_cur(T)); chk('TT/5 = 0', TT[5], {})
    for nm, A, h in (('J', J, 1), ('G+', Gp2, sp.Rational(3, 2)), ('G-', Gm2, sp.Rational(3, 2))):
        TA = poles(T, A)
        chk(f'T{nm}/3 = 0 (primary)', TA[3], {}); chk(f'T{nm}/4 = 0', TA[4], {})
        chk(f'T{nm}/2 = {h} {nm}', TA[2], scale(A, h)); chk(f'T{nm}/1 = d{nm}', TA[1], d_cur(A))
    JJ = poles(J, J)
    log(f"  J J double pole = {sp.factor(JJ[2].get((), 0))}  (c/3 = {sp.factor(c/3)})")
    chk('JJ/2 = c/3', JJ[2], scale(one, c / 3)); chk('JJ/1 = 0', JJ[1], {}); chk('JJ/3 = 0', JJ[3], {})
    JGp = poles(J, Gp2); JGm = poles(J, Gm2)
    sgn = None
    for sg in (1, -1):
        if all(sp.simplify(JGp[1].get(m, 0) - sg * Gp2.get(m, 0)) == 0 for m in set(JGp[1]) | set(Gp2)):
            sgn = sg
    log(f"  J G+ simple pole = {sgn} G+ (sign convention of J-charge)")
    chk('JG+/1 = +-G+', JGp[1], scale(Gp2, sgn if sgn else 1)); chk('JG+/2 = 0', JGp[2], {})
    chk('JG-/1 = -+G-', JGm[1], scale(Gm2, -sgn if sgn else -1)); chk('JG-/2 = 0', JGm[2], {})
    PP = poles(Gp2, Gp2); MM = poles(Gm2, Gm2); PM = poles(Gp2, Gm2)
    for n in (1, 2, 3):
        chk(f'G+G+/{n} = 0', PP[n], {}); chk(f'G-G-/{n} = 0', MM[n], {})
    # Gp2 Gm2 = 2 G+ G-: 2[(2c/3)/z^3 + 2J/z^2 + (2T + dJ)/z]
    chk('G+G-/3 = 2c/3 (x2)', PM[3], scale(one, 4 * c / 3)); chk('G+G-/2 = 2J (x2)', PM[2], scale(J, 4))
    chk('G+G-/1 = 2T + dJ (x2)', PM[1], add((4, T), (2, d_cur(J))))
    log(f"V1 result: c = {c}; all N = 2 OPEs (complex basis T, G+, G-, J; KL (9) form): {'ALL OK' if allok else 'FAILURES'}  ({time.time()-t0:.1f}s)")
    return dict(ok=bool(allok), c=str(c), checks=res, gprime={str(k): str(v) for k, v in gp.items()},
                J={str(k): str(v) for k, v in J.items()}, J_charge_sign_of_Gplus=sgn)


# ------------------------------------------------------------------ vertices

def kl_vertices(sv):
    """KL S1..S4 as SUSY vertices (b.psi) e^{b.X}: S1 = m_1, S2 = m_2, S3 = -m_1, S4 = -m_2 (labels as KL (16), (19)
    up to the overall normalisation of the prefactor)."""
    h = F(1, 2)
    m = {1: (h, h), 2: (h, -h), 3: (-h, -h), 4: (-h, h)}
    return {k: susy_vertex(v, eta_of(sv), label=f'S{k}') for k, v in m.items()}


def w0_vertex(sv, use='G'):
    """Zeitlin's third vertex W_0 = int dtheta (D Phi_1^Z - D Phi_2^Z) e^{Phi_1^Z + Phi_2^Z} in our frame:
    Phi_1^Z = -(Phi_0 + Phi_1)/2, Phi_2^Z = -(Phi_0 - Phi_1)/2, so W_0 = -int dtheta (D Phi_1) e^{-Phi_0}.
    Bottom component B = psi_1 e^{-X_0}; top component = G_{-1/2} B = Res(G, B) (lead's instruction, item 348).
    Returns (bottom Vertex, top Vertex)."""
    C = currents(sv)
    Gsel = C['G'] if use == 'G' else C['gp']
    eta = eta_of(sv)
    b = (F(-1), F(0))
    bottom = Vertex(b, {(PSI(0, 1),): F(1)}, label='W0_bottom')
    top_W = vertex_residue(Gsel, bottom, eta)
    top = Vertex(b, top_W, label=f'W0_top[{use}]')
    return bottom, top


def task_weights(log):
    eta = eta_of(s)
    C = currents(s)
    T = C['T']
    out = {}
    for k, V in kl_vertices(s).items():
        Wsp = {m: sp.Rational(v.numerator, v.denominator) for m, v in V.W.items()}
        Vs = Vertex(V.b, Wsp, label=V.label)
        hp = vertex_residue(T, Vs, eta, shift=1)
        m0 = next(iter(Wsp))
        hv = sp.simplify(hp.get(m0, 0) / Wsp[m0])
        prop = all(sp.simplify(hp.get(m, 0) - hv * Wsp.get(m, 0)) == 0 for m in set(hp) | set(Wsp))
        prim = all(not vertex_residue(T, Vs, eta, shift=sh) for sh in (2, 3))
        out[V.label] = dict(h=str(hv), proportional=prop, primary=prim)
        log(f"  weight {V.label}: h = {hv}; double pole proportional: {prop}; primary: {prim}")
    # W_0: exponential part and top component
    bottom, top = w0_vertex(s)
    for V in (boson_vertex((F(-1), F(0)), label='e^{-X_0}'), bottom, top):
        hp = vertex_residue(T, V, eta, shift=1)
        Wd = {m: sp.sympify(v) for m, v in V.W.items()}
        m0 = next(iter(Wd))
        hv = sp.simplify(hp.get(m0, 0) / Wd[m0]) if hp.get(m0, 0) != 0 else None
        prop = hv is not None and all(sp.simplify(hp.get(m, 0) - hv * Wd.get(m, 0)) == 0 for m in set(hp) | set(Wd))
        higher = [vertex_residue(T, V, eta, shift=sh) for sh in (2, 3)]
        out[V.label] = dict(h=str(hv), proportional=prop, primary=all(not h_ for h_ in higher),
                            prefactor={str(k): str(v) for k, v in V.W.items()})
        log(f"  weight {V.label}: h = {hv}; proportional: {prop}; higher poles vanish: {all(not h_ for h_ in higher)}; prefactor {V.W}")
    return out


def task_vertices(log, svals=(F(2, 3), F(10, 9), F(8, 7), F(8, 9), F(3), F(5, 7))):
    """Commutation of T, G, G', J with S1..S4 and with W_0 (top, built from G and from G') at rational s."""
    out = {}
    for sv in svals:
        C = currents(sv)
        eta = C['eta']
        S = kl_vertices(sv)
        rows = {}
        for nm in ('T', 'G', 'gp', 'J'):
            rows[nm] = {f'S{k}': screens(C[nm], S[k], eta) for k in (1, 2, 3, 4)}
        for use in ('G', 'gp'):
            bottom, top = w0_vertex(sv, use)
            for nm in ('T', 'G', 'gp', 'J'):
                rows[nm][f'W0[{use}]'] = screens(C[nm], top, eta)
            rows.setdefault('bottom', {})[f'W0b[{use}]'] = {nm: screens(C[nm], bottom, eta) for nm in ('G', 'gp', 'J')}
        out[str(sv)] = rows
        log(f"  s = {sv}: " + ' | '.join(f"{nm}: " + ' '.join(f"{k}={'Y' if v else 'n'}" for k, v in r.items())
                                        for nm, r in rows.items() if nm != 'bottom'))
        log(f"      bottom psi_1 e^(-X_0): {rows['bottom']}")
    return out


if __name__ == '__main__':
    task = sys.argv[1]
    lines = []

    def log(msg):
        print(msg, flush=True)
        lines.append(msg)
    t0 = time.time()
    if task == 'v1':
        r = task_v1(log)
    elif task == 'weights':
        r = task_weights(log)
    elif task == 'vertices':
        r = task_vertices(log)
    else:
        raise SystemExit('task?')
    json.dump(r, open(f'n2frame_{task}.json', 'w'), indent=1, default=str)
    print(f'done {task} in {time.time()-t0:.1f}s')
