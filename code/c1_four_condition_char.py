"""The four-condition CHARACTERIZATION of the enhanced c=1 tower.

Upgrade prompted by adversarial review: the four momenta
{+-1/sqrt2, +-2 sqrt2} are, at c=1, exactly the Kac momenta
alpha_{r,s} = (s-r)/sqrt2 of the degenerate fields Phi_{1,2},
Phi_{2,1} (h = 1/4) and Phi_{1,5}, Phi_{5,1} (h = 4).  This
certificate proves the sharp converse of the fingerprint table:
the joint kernel of ALL FOUR first-order conservation conditions
has charge dimension exactly ONE at every odd spin 1..9 and ZERO
at every even spin through the range -- so, through the computed
range, the common conserved structure of the Phi_{12}/Phi_{21} and
Phi_{15}/Phi_{51} perturbations at c=1 IS the enhanced tower
(together with oint T), with multiplicity one.  The four
conditions do not merely admit the tower; they characterize it.

Honest scoping (corrected on review): ALL FOUR labels belong to
the single known A2(2) family of integrable directions --
Phi_12/Phi_21 classically (Zamolodchikov, Smirnov) and Phi_15 as
an established integrable perturbation related to A2(2) and the
Izergin-Korepin model (Kausch-Takacs-Watts NPB489 (1997) 557;
Berkovich-McCoy-Pearce NPB519 (1998) 597).  Those are c<1
M(p,p') theorems; at c=1 the Kac labels are formal degenerate
weights, so the classical results are the antecedent, not
transferable by relabelling.  Relevance (h=4 at c=1) bears on
flows, not on the first-order condition.
Restriction of scalars over Q(sqrt2), sparse integer elimination;
dims are exact.
"""
import sys
from fractions import Fraction as F
import sympy as sp
from ope_engine import d_mono, gen_basis
from screening_fit import residue_screen, twisted_d
from sparse_linalg import SparsePivots, lcm

SQ2 = sp.sqrt(2)
NB = 2

def split2(expr):
    e = sp.expand(expr)
    c1, c2 = e.coeff(SQ2, 0), e.coeff(SQ2, 1)
    if sp.simplify(e - c1 - c2 * SQ2) != 0:
        raise ValueError(f"split reconstruction failed: {e}")
    r1, r2 = sp.Rational(c1), sp.Rational(c2)
    return (F(int(r1.p), int(r1.q)), F(int(r2.p), int(r2.q)))

def mul2(pair, bi):
    c1, c2 = pair
    return (c1, c2) if bi == 0 else (2 * c2, c1)

def int_vec(fr):
    L = 1
    for v in fr.values():
        L = lcm(L, v.denominator)
    return {i: int(v * L) for i, v in fr.items() if v}, L

def joint(w, moms):
    basisw = gen_basis(w, 1, z2=False)
    tgt = gen_basis(w - 1, 1, z2=False)
    tidx = {m: i for i, m in enumerate(tgt)}
    lower = [()] if w == 2 else gen_basis(w - 2, 1, z2=False)
    nt = len(tgt)
    R = [[{tidx[m]: split2(c) for m, c in residue_screen({b: F(1)}, al).items()} for b in basisw] for al in moms]
    T = [[{tidx[m]: split2(c) for m, c in twisted_d(u, al, 1).items()} for u in lower] for al in moms]
    cols, scales = [], []
    nx = len(basisw) * NB
    for j in range(len(basisw)):
        for bi in range(NB):
            col = {}
            for s in range(len(moms)):
                base = s * nt * NB
                for r, e in R[s][j].items():
                    for beta, cv in enumerate(mul2(e, bi)):
                        if cv:
                            col[base + r*NB + beta] = col.get(base + r*NB + beta, F(0)) + cv
            iv, L = int_vec(col); cols.append(iv); scales.append(L)
    for s in range(len(moms)):
        base = s * nt * NB
        for u in range(len(lower)):
            for bi in range(NB):
                col = {}
                for r, e in T[s][u].items():
                    for beta, cv in enumerate(mul2(e, bi)):
                        if cv:
                            col[base + r*NB + beta] = col.get(base + r*NB + beta, F(0)) + cv
                iv, L = int_vec(col); cols.append(iv); scales.append(L)
    nc = len(cols)
    piv = SparsePivots(); kers = []
    for jc, col in enumerate(cols):
        aug = {i + nc: v for i, v in col.items()}; aug[jc] = 1
        r = piv.reduce(aug)
        if r and max(r) < nc: kers.append(r)
        elif r: piv.rows[max(r)] = r
    vecs = []
    for kv in kers:
        xv = {c: co * scales[c] for c, co in kv.items() if c < nx and co * scales[c]}
        if xv: vecs.append(xv)
    # mod-d
    basis = gen_basis(w, 1, z2=False); idx = {m: i for i, m in enumerate(basis)}
    dp = SparsePivots()
    for bm in gen_basis(w - 1, 1, z2=False):
        col = {}
        for m, c in d_mono(bm).items(): col[idx[m]] = col.get(idx[m], 0) + c
        col = {i: c for i, c in col.items() if c}
        if col:
            for beta in range(NB): dp.insert({i*NB+beta: c for i, c in col.items()})
    qd = sum(1 for v in vecs if v and dp.insert(dict(v)))
    return qd

failures = []


def check(name, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {name}", flush=True)
    if not ok:
        failures.append(name)


MOMS4 = [(SQ2/2,), (-SQ2/2,), (2*SQ2,), (-2*SQ2,)]
for (r, s_), al_want, h_want in (((1, 2), SQ2/2, sp.Rational(1, 4)),
                                 ((2, 1), -SQ2/2, sp.Rational(1, 4)),
                                 ((1, 5), 2*SQ2, sp.Integer(4)),
                                 ((5, 1), -2*SQ2, sp.Integer(4))):
    al = sp.Rational(s_ - r) / SQ2
    check(f"Kac dictionary: alpha_{{{r},{s_}}} = (s-r)/sqrt2 = "
          f"{al_want} with h = {h_want}",
          sp.simplify(al - al_want) == 0
          and sp.simplify(al**2/2 - h_want) == 0)

EXPECT = {2: 2, 3: 0, 4: 2, 5: 0, 6: 2, 7: 0, 8: 2, 9: 0, 10: 2}
for w in range(2, 11):
    qd = joint(w, MOMS4)
    check(f"weight {w} (spin {w-1}): joint-kernel charge Q-dim = "
          f"{EXPECT[w]} (K-dim {EXPECT[w]//2})", qd == EXPECT[w])

print()
if failures:
    print(f"FAILURES: {failures}")
    print("FOUR-CONDITION CHARACTERIZATION FAILED")
    sys.exit(1)
print("FOUR-CONDITION CHARACTERIZATION CERTIFIED: through weights "
      "2..10 the joint kernel of the Phi_12/Phi_21/Phi_15/Phi_51 "
      "first-order conservation conditions at c=1 has charge "
      "multiplicity ONE at every odd spin and ZERO at every even "
      "spin -- the enhanced tower (with the spin-1 charge "
      "oint T) IS the common "
      "conserved structure, exactly.")
print("FOUR-CONDITION CHARACTERIZATION COMPLETE")
