"""classical.py -- B2 V2: the classical limit of Towers Z and K against the N = 2 KdV family SKdV_a.

(1) Classical commutants: the screening conditions with ONE Wick contraction (gm.joint_kernel(classical=True)) at the
    classical metric eta-hat = (1, -1) (the quantum metric is s eta-hat; a k-contraction term carries s^k).  The classical
    currents are n2frame.currents(1): T^ = 1/2(dX_0)^2 - 1/2(dX_1)^2 - 1/2 psi_0 dpsi_0 + 1/2 psi_1 dpsi_1 - d^2 X_0,
    J^ = psi_0 psi_1 - 2 dX_1, G^, g'^ (= s T, s J, ... of the quantum ones; Zeitlin's (5) up to conventions).
(2) Quantum -> classical check at w2 = 8: the quantum kernel density (normalised: coefficient of (dX_0)^4 = 1) at
    s = 1/j, j = 2..26, each coefficient reconstructed as a rational function of s (exact linear algebra, verified at
    5 further points) and evaluated at s = 0; must equal the classical kernel density.
(3) Flow: Phi_t = {oint D_cl, Phi} with the single-contraction bracket, for the N = 2 superfield
    Phi = alpha J^ + theta_1 (b1 G+ + b2 G-) + theta_2 (g1 G+ + g2 G-) + theta_1 theta_2 (d T^ + e dJ^),
    whose superfield structure (b, g, d, e) is first fixed by requiring {oint G^, .} and {oint g'^, .} to act as k_1 D_1
    and k_2 D_2 (both assignments tried); then compared with Mao-Liu 1705.03997 (1) (= Laberge-Mathieu):
    phi_t = -phi_xxx + 3(phi D_1 D_2 phi)_x + (a-1)/2 (D_1 D_2 phi^2)_x + 3a phi^2 phi_x,  D_i = d/dtheta_i + theta_i d_x,
    phi = v + theta_1 xi_1 + theta_2 xi_2 + theta_1 theta_2 u (components defined by this expansion), with the time scale
    tau and a as unknowns.  Output: a for each tower (or "no a").
Usage: classical.py <tower Z|K> [w2max]
"""
import sys
sys.dont_write_bytecode = True
import json
import time
from fractions import Fraction as F

import sympy as sp

import gm
from gm import X, PSI, add, scale, d_cur, clean
from super2_screening import Vertex, sector_candidates, susy_vertex
import n2frame as N2
from towers import cur_from_json

ETA_CL = (F(1), F(-1))


def classical_vertices(tower):
    S = N2.kl_vertices(F(1))
    C = N2.currents(F(1))
    if tower == 'K':
        return [S[1], S[2], S[3], S[4]]
    bottom = Vertex((F(-1), F(0)), {(PSI(0, 1),): F(1)}, label='W0_bottom')
    top_W = gm.vertex_residue1(C['G'], bottom, ETA_CL)
    return [S[1], S[2], Vertex((F(-1), F(0)), top_W, label='W0_cl')]


# ------------------------------------------------------------------ classical products and superfields

def cmul(A, B):
    """classical (concatenation) product of letter polynomials, Grassmann signs by canonicalize."""
    from super_engine import canonicalize
    out = {}
    for ma, ca in A.items():
        for mb, cb in B.items():
            cm, sg = canonicalize(ma + mb)
            if cm is None:
                continue
            out[cm] = out.get(cm, 0) + sg * ca * cb
    return clean(out)


def parity(cur):
    ps = {sum(1 for l in m if l[0] == 1) % 2 for m in cur}
    assert len(ps) <= 1, ps
    return ps.pop() if ps else 0


TH = [(), (1,), (2,), (1, 2)]


def th_mul(a, b):
    w = list(a) + list(b)
    if len(set(w)) != len(w):
        return 0, None
    sign = 1
    w2 = w[:]
    for i in range(len(w2)):
        for j in range(len(w2) - 1 - i):
            if w2[j] > w2[j + 1]:
                w2[j], w2[j + 1] = w2[j + 1], w2[j]
                sign = -sign
    return sign, tuple(w2)


def sf_clean(sf):
    return {th: clean(c) for th, c in sf.items() if clean(c)}


def sf_add(*terms):
    out = {}
    for coef, sf in terms:
        for th, c in sf.items():
            out[th] = add((1, out.get(th, {})), (coef, c))
    return sf_clean(out)


def sf_mul(A, B):
    """product of superfields (components are letter polynomials; theta's anticommute with odd components)."""
    out = {}
    for tha, ca in A.items():
        pa = parity(ca)
        for thb, cb in B.items():
            sg, th = th_mul(tha, thb)
            if not sg:
                continue
            sg *= (-1) ** (pa * len(thb))
            out[th] = add((1, out.get(th, {})), (sg, cmul(ca, cb)))
    return sf_clean(out)


def sf_D(A, i):
    """D_i = d/dtheta_i + theta_i d_x (acting from the left; components sit to the right of the theta's)."""
    out = {}
    for th, c in A.items():
        if i in th:
            pos = th.index(i)
            sg = (-1) ** pos
            new = tuple(x for x in th if x != i)
            out[new] = add((1, out.get(new, {})), (sg, c))
        else:
            sg = (-1) ** sum(1 for x in th if x < i)
            new = tuple(sorted(th + (i,)))
            out[new] = add((1, out.get(new, {})), (sg, d_cur(c)))
    return sf_clean(out)


def sf_Q(A, i):
    """Q_i = d/dtheta_i - theta_i d_x (the superspace supercharge; [Q_i, D_j} = 0)."""
    out = {}
    for th, c in A.items():
        if i in th:
            pos = th.index(i)
            sg = (-1) ** pos
            new = tuple(x for x in th if x != i)
            out[new] = add((1, out.get(new, {})), (sg, c))
        else:
            sg = -(-1) ** sum(1 for x in th if x < i)
            new = tuple(sorted(th + (i,)))
            out[new] = add((1, out.get(new, {})), (sg, d_cur(c)))
    return sf_clean(out)


def sf_dx(A):
    return sf_clean({th: d_cur(c) for th, c in A.items()})


def sf_scale(A, k):
    return sf_clean({th: scale(c, k) for th, c in A.items()})


def bracket(H, Fcur):
    """{oint H, F(w)} classical = single-contraction residue (metric eta-hat)."""
    return gm.residue1(H, Fcur, ETA_CL)


def sf_bracket(H, A):
    """{oint H, Phi} for a superfield Phi = sum theta^th A_th.  For an ODD charge H the graded bracket is taken with
    the Grassmann parameter pulled to the left, [eps oint H, theta^th A] = eps (-1)^|th| theta^th [oint H, A}, so the
    components behind an odd number of theta's pick up a sign; an even charge gives no sign."""
    hp = parity(H)
    return sf_clean({th: scale(bracket(H, c), (-1) ** (hp * len(th))) for th, c in A.items()})


# ------------------------------------------------------------------ rational reconstruction in s

def reconstruct(samples, extra, maxdeg=10):
    """samples: list of (s, value) Fractions.  Find num/den of minimal total degree reproducing all samples, verify on
    extra; return (num coeffs, den coeffs) or None."""
    for dn in range(0, maxdeg + 1):
        for dd in range(0, maxdeg + 1 - dn):
            nun = dn + 1 + dd           # unknowns: n_0..n_dn, d_1..d_dd (d_0 = 1)
            if len(samples) < nun + 1:
                continue
            rows = []
            rhs = []
            for sv, fv in samples:
                rows.append([sv ** i for i in range(dn + 1)] + [-fv * sv ** i for i in range(1, dd + 1)])
                rhs.append(fv)
            M = sp.Matrix(rows)
            r = sp.Matrix(rhs)
            try:
                sol, params = M.gauss_jordan_solve(r)
            except ValueError:
                continue
            if params:
                continue
            num = [sol[i] for i in range(dn + 1)]
            den = [sp.Integer(1)] + [sol[dn + 1 + i] for i in range(dd)]
            ok = True
            for sv, fv in extra:
                n = sum(num[i] * sv ** i for i in range(dn + 1))
                d = sum(den[i] * sv ** i for i in range(dd + 1))
                if d == 0 or sp.Rational(n) / d != fv:
                    ok = False
                    break
            if ok:
                return num, den
    return None


def main():
    tower = sys.argv[1]
    w2max = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    log = []

    def P(msg):
        print(msg, flush=True)
        log.append(msg)
    t0 = time.time()
    V = classical_vertices(tower)
    P(f"classical tower {tower}: vertices " + ', '.join(f"{v.label}: b={v.b}, W={v.W}" for v in V))
    C = N2.currents(F(1))
    dims = {}
    kers = {}
    for w2 in range(2, w2max + 1):
        ker = gm.joint_kernel(sector_candidates(w2, 'all'), V, ETA_CL, classical=True)
        dims[w2] = len(ker)
        kers[w2] = ker
        P(f"  classical commutant w2 = {w2}: dim {len(ker)}")
    out = dict(tower=tower, classical_dims=dims)
    # which low charges: J^, G^, g'^, T^ in the classical kernel?
    for nm in ('J', 'G', 'gp', 'T'):
        cur = C[nm]
        w2 = gm.mono_weight2(next(iter(cur)))
        P(f"  {nm}^ in classical commutant (mod d): {gm.in_span_mod_d(cur, kers.get(w2, []))}")
    assert dims.get(8) == 1, dims
    Dcl = kers[8][0]
    top = (X(1, 0),) * 4
    Dcl = scale(Dcl, F(1) / Dcl[top])
    P(f"classical weight-4 density (normalised (dX_0)^4 = 1): {len(Dcl)} monomials")
    out['D_cl'] = [[list(map(list, m)), str(c)] for m, c in sorted(Dcl.items())]

    # (2) quantum -> classical limit at w2 = 8 by rational reconstruction
    P("quantum kernels at s = 1/j, j = 2..31 (w2 = 8) ...")
    quantum = {}
    for j in range(2, 32):
        sv = F(1, j)
        Vq = N2.kl_vertices(sv)
        if tower == 'K':
            Vl = [Vq[1], Vq[2], Vq[3], Vq[4]]
        else:
            Vl = [Vq[1], Vq[2], N2.w0_vertex(sv, 'G')[1]]
        ker = gm.joint_kernel(sector_candidates(8, 'all'), Vl, N2.eta_of(sv))
        assert len(ker) == 1, (sv, len(ker))
        k = scale(ker[0], F(1) / ker[0][top])
        quantum[sv] = k
    monos = sorted(set(m for k in quantum.values() for m in k) | set(Dcl))
    limit = {}
    allok = True
    for m in monos:
        pts = [(sp.Rational(sv), sp.Rational(k.get(m, 0))) for sv, k in quantum.items()]
        rec = reconstruct(pts[:25], pts[25:])
        if rec is None:
            P(f"  reconstruction FAILED for {m}")
            allok = False
            continue
        num, den = rec
        s0 = num[0] / den[0]
        limit[m] = s0
        if sp.Rational(s0) != sp.Rational(Dcl.get(m, 0)):
            P(f"  s -> 0 limit MISMATCH at {m}: limit {s0} vs classical {Dcl.get(m, 0)}; num {num} den {den}")
            allok = False
    P(f"quantum (s -> 0, reconstructed from 25 points, verified on 5) == classical commutant density: {allok}")
    out['limit_equals_classical'] = allok

    # (3) superfield structure and the flow
    T, G, gp, J = C['T'], C['G'], C['gp'], C['J']
    Gp2, Gm2 = C['Gp2'], C['Gm2']
    b1, b2, g1, g2, dd, ee, k1, k2 = sp.symbols('b1 b2 g1 g2 d e k1 k2')
    Jsp = {m: sp.Rational(c) for m, c in J.items()}
    Tsp = {m: sp.Rational(c) for m, c in T.items()}
    Gpsp = {m: sp.Rational(c) for m, c in Gp2.items()}
    Gmsp = {m: sp.Rational(c) for m, c in Gm2.items()}
    Phi = {(): Jsp, (1,): add((b1, Gpsp), (b2, Gmsp)), (2,): add((g1, Gpsp), (g2, Gmsp)),
           (1, 2): add((dd, Tsp), (ee, d_cur(Jsp)))}
    Phi = sf_clean(Phi)
    results = {}
    for gen, assign in [(g, a_) for g in ('Q', 'D') for a_ in (('G', 'gp'), ('gp', 'G'))]:
        Q1 = {m: sp.Rational(c) for m, c in C[assign[0]].items()}
        Q2 = {m: sp.Rational(c) for m, c in C[assign[1]].items()}
        eqs = []
        for Q, i, kk in ((Q1, 1, k1), (Q2, 2, k2)):
            lhs = sf_bracket(Q, Phi)
            rhs = sf_scale((sf_Q if gen == 'Q' else sf_D)(Phi, i), kk)
            diff = sf_add((1, lhs), (-1, rhs))
            for th, c in diff.items():
                eqs += [sp.expand(v) for v in c.values()]
        eqs = [e for e in eqs if e != 0]
        sols = sp.solve(eqs, [b1, b2, g1, g2, dd, ee, k1, k2], dict=True)
        P(f"superfield structure [{gen}_i] with ({assign[0]} -> {gen}_1, {assign[1]} -> {gen}_2): {len(sols)} solution(s): {sols}")
        for sol in sols:
            if any(sol.get(x, 1) == 0 for x in (dd, k1, k2)):
                continue
            Phi_s = sf_clean({th: {m: sp.expand(sp.sympify(v).subs(sol)) for m, v in c.items()} for th, c in Phi.items()})
            # flow of Phi under H = oint D_cl
            Dsp = {m: sp.Rational(c) for m, c in Dcl.items()}
            lhs = sf_bracket(Dsp, Phi_s)
            a, tau, al = sp.symbols('a tau alpha')
            # Mao-Liu's phi = alpha * Phi_s (the equation is not homogeneous in phi: alpha is an unknown)
            D12 = sf_D(sf_D(Phi_s, 2), 1)
            Phi2 = sf_mul(Phi_s, Phi_s)
            D12_2 = sf_D(sf_D(Phi2, 2), 1)
            rhs = sf_add((-1, sf_dx(sf_dx(sf_dx(Phi_s)))), (3 * al, sf_dx(sf_mul(Phi_s, D12))),
                         (al * (a - 1) / 2, sf_dx(D12_2)), (3 * a * al ** 2, sf_mul(Phi2, sf_dx(Phi_s))))
            diff = sf_add((tau, lhs), (-1, rhs))
            eqs2 = [sp.expand(v) for c in diff.values() for v in c.values()]
            eqs2 = [e for e in eqs2 if e != 0]
            sol2 = sp.solve(eqs2, [a, tau, al], dict=True)
            P(f"   flow fit (Mao-Liu (1), phi = alpha Phi_s) with {sol}: {sol2 if sol2 else 'NO (a, tau, alpha) solves all component equations'}")
            if not sol2:
                # diagnose: least-squares-free: show how many equations and a sample
                P(f"      {len(eqs2)} equations; sample: {eqs2[:3]}")
            results[f"{gen}:{assign}:{sol}"] = str(sol2)
    out['flow_fits'] = results
    json.dump(out, open(f'classical_{tower}.json', 'w'), indent=1, default=str)
    P(f"done in {time.time()-t0:.1f}s")


if __name__ == '__main__':
    main()
