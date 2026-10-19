"""Validation of super_cyl_vev.py (SUSY Phase 4, Track T1).  All exact.

 (i)   T_l -> L_0 - c/24 in NS and R, symbolic in P and l, with c(l) read
       off the ENGINE (fourth-order pole of T_l T_l, full contractions with
       super_engine's pairings and Koszul signs) and the Fock dimensions from
       the mode algebra: plane L_0 = a_0^2/2 - l a_0 + (oscillators)
       [z^2 d^2X has zero mode -a_0], R adds the Majorana zero-point 1/16
       [L_0^f = sum_{n>0} n psi_{-n} psi_n + 1/16]; evaluator momentum
       P = p - l.  Expected: NS  P^2/2 - 1/16,  R  P^2/2.
 (ii)  R sector: G_0^2 = L_0 - c/24 on the R ground states, G_0 evaluated
       as an element of the zero-mode Clifford algebra (psi_0^2 = 1/2).
 (iii) negative controls, each MUST FAIL the check it targets:
       NS<->R fermion kernels swapped; a perturbed density coefficient;
       a perturbed kernel derivative; the R zero mode dropped.
 (iv)  purely bosonic monomials (w2 <= 16): identical to cyl_vev.py
       (exact Fractions, N = 1, Q = 0), in both sectors.
 (v)   total derivatives have zero eigenvalue: vev(d b) = 0 for every
       canonical monomial b with w2(b) <= 12, both sectors, both parities.
Exit status 0 iff every check behaves as required."""
import sys
from fractions import Fraction as F
from itertools import permutations

import sympy as sp

import cyl_vev
import super_cyl_vev as scv
from super_engine import _koszul, _pairing, d_mono, gen_basis_super
from super_lambda import G_l, T_l

fails = []


def require(cond, msg):
    print(("PASS  " if cond else "FAIL  ") + msg, flush=True)
    if not cond:
        fails.append(msg)


def must_fail(cond_ok, msg):
    """negative control: the targeted check must NOT hold"""
    print(("PASS  " if not cond_ok else "FAIL  ") + "[neg. control fails as required] " + msg, flush=True)
    if cond_ok:
        fails.append("negative control did not fail: " + msg)


def central_term(A, B, order):
    """c-number coefficient of (z-w)^-order in A(z)B(w) (full contractions)."""
    tot = 0
    for ma, ca in A.items():
        for mb, cb in B.items():
            if len(ma) != len(mb):
                continue
            n = len(ma)
            for qsel in permutations(range(n)):
                num, D = 1, 0
                for a, b in zip(range(n), qsel):
                    pr = _pairing(ma[a], mb[b])
                    if pr is None:
                        num = 0
                        break
                    num *= pr[0]
                    D += pr[1]
                if num and D == order:
                    tot += ca * cb * num * _koszul(ma, mb, tuple(range(n)), qsel)
    return sp.expand(tot)


l, P = sp.symbols('l P')


def check_T(sector, Tdens=None, **kw):
    Tdens = Tdens if Tdens is not None else T_l(l)
    got = scv.scalar_vev(Tdens, sector, [P], **kw)
    c = 2 * central_term(T_l(l), T_l(l), 4)
    p = P + l
    Delta = p ** 2 / 2 - l * p + (sp.Rational(1, 16) if sector == 'R' else 0)
    return sp.expand(got - (Delta - c / 24)) == 0, got, c


def check_G0sq(**kw):
    g = scv.vev_density(G_l(l), 'R', [P], **kw)
    odd = g.get((0,), sp.Integer(0))
    even = g.get((), sp.Integer(0))
    G0sq = sp.expand(odd ** 2 * sp.Rational(1, 2))       # (b psi_0)^2 = b^2/2
    t = scv.scalar_vev(T_l(l), 'R', [P], **kw)
    return even == 0 and sp.expand(G0sq - t) == 0, odd, t


if __name__ == "__main__":
    # (i)
    for sec in ('NS', 'R'):
        ok, got, c = check_T(sec)
        require(ok, f"(i) {sec}: <T_l> = {got} = Delta_{sec}(p) - c/24 with engine c(l) = {c}, P = p - l")
    require(sp.expand(2 * central_term(T_l(l), T_l(l), 4) - (sp.Rational(3, 2) - 12 * l ** 2)) == 0,
            "(i) engine central charge c(l) = 3/2 - 12 l^2 (= Tanabe (4.16) at alpha_0 = -2l)")
    # (ii)
    ok, odd, t = check_G0sq()
    require(ok, f"(ii) R: G_0 = ({odd}) psi_0 on the ground space, G_0^2 = {sp.expand(odd**2/2)} = L_0 - c/24 = {t}")
    # (iii) negative controls
    swapped = {'NS': ('B', 'R'), 'R': ('B', 'NS')}
    must_fail(check_T('NS', kernels=swapped)[0], "(iii) NS check with the R kernel")
    must_fail(check_T('R', kernels=swapped)[0], "(iii) R check with the NS kernel")
    pert = dict(T_l(l))
    key = [m for m in pert if all(x[0] == 1 for x in m)][0]
    pert[key] = pert[key] + sp.Rational(1, 100)
    must_fail(check_T('NS', Tdens=pert)[0], "(iii) NS check with the fermion coefficient of T_l perturbed by 1/100")
    must_fail(check_T('R', Tdens=pert)[0], "(iii) R check with the fermion coefficient of T_l perturbed by 1/100")
    # perturbed kernel derivative: K_NS'(0) -> K_NS'(0) + 1/1000 via a custom table
    orig = scv.kernel_table('NS')
    scv.kernel_table.cache_clear()
    scv._fermion_part.cache_clear()
    _kt = scv.kernel_table

    def bad_table(kind, n=scv._NMAX, _o=orig):
        if kind == 'NSbad':
            t = list(_o)
            t[1] += F(1, 1000)
            return tuple(t)
        return _kt(kind, n)
    scv.kernel_table = bad_table
    must_fail(check_T('NS', kernels={'NS': ('B', 'NSbad'), 'R': ('B', 'R')})[0],
              "(iii) NS check with K_NS'(0) perturbed by 1/1000")
    scv.kernel_table = _kt
    scv._fermion_part.cache_clear()
    must_fail(check_G0sq(keep_zero_mode=False)[0], "(iii) R G_0^2 check with the zero mode psi_0 dropped")
    # (iv) bosonic sector identical to cyl_vev.py
    nb, bad = 0, []
    for w2 in range(2, 17, 2):
        for mono in gen_basis_super(w2, 1):
            if any(x[0] == 1 for x in mono):
                continue
            xs = tuple(x[1] for x in mono)
            for sec in ('NS', 'R'):
                r = scv.vev_monomial(mono, sec)
                # cyl_vev.py is the O(N) evaluator: momentum enters only as P^2 and
                # odd numbers of unpaired zero modes are discarded by construction
                # (its `nx % 2` rule).  Compare the even-in-P part, which is the
                # whole answer for monomials with an even number of letters.
                d = {e: v for e, v in r.get((), {}).items() if e[0] % 2 == 0}
                if len(mono) % 2 == 0 and len(d) != len(r.get((), {})):
                    bad.append((mono, sec, 'odd power in an even monomial'))
                for P2 in (F(0), F(1), F(2), F(-3), F(1, 3), F(7, 5), F(-11, 4)):
                    ours = sum((v * P2 ** (e[0] // 2) for e, v in d.items()), F(0))
                    theirs = cyl_vev.vev_monomial(xs, (), F(1), P2, F(0))
                    if ours != theirs:
                        bad.append((mono, sec, P2, ours, theirs))
            nb += 1
    require(not bad, f"(iv) {nb} bosonic monomials (w2 <= 16), NS and R, at 7 momenta: even-in-P part identical to cyl_vev.py exactly"
            + (f"; mismatches {bad[:3]}" if bad else ""))
    # (v) total derivatives
    nt, bad = 0, []
    for w2 in range(1, 15):
        for b in gen_basis_super(w2 - 2, 1) if w2 >= 3 else []:
            db = d_mono(b)
            for sec in ('NS', 'R'):
                if scv.vev_density(db, sec, [P]):
                    bad.append((b, sec))
            nt += 1
    require(not bad, f"(v) vev(d b) = 0 for all {nt} monomials b (w2(b) <= 12), NS and R" + (f"; failures {bad[:3]}" if bad else ""))
    print()
    if fails:
        print(f"SUPER CYL VEV VALIDATION: {len(fails)} FAILURE(S)")
        sys.exit(1)
    print("SUPER CYL VEV VALIDATED")
