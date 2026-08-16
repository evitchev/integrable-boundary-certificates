"""V1'-suite for the N=1 super-OPE engine (super_engine.py).

Mirrors validate.py's discipline: exact checks, a decisive canary, a
perturbed-coefficient negative control, and a cross-check against the
trusted bosonic engine on the purely bosonic subalgebra.

The engine is residue-only, so central terms ((z-w)^-4 in TT,
(z-w)^-3 in GG) are invisible here by design; they do not enter charge
brackets.  NS/R sector choices never enter these computations.
"""
import sys
from fractions import Fraction as F

import ope_engine
from super_engine import (SuperSector, canonicalize, cur_add, d_cur,
                          d_mono, gen_basis_super, mono_weight2,
                          normal_prod, residue_cur, residue_mono)

DX = (0, 1, 0)
D2X = (0, 2, 0)
PSI = (1, 0, 0)
DPSI = (1, 1, 0)

# T = 1/2 (dX)^2 - 1/2 psi dpsi   (c = 3/2), G = psi dX
T = {(DX, DX): F(1, 2), (PSI, DPSI): F(-1, 2)}
G = {(DX, PSI): F(1)}

failures = []


def check(name, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {name}")
    if not ok:
        failures.append(name)


def res_cur_mono(cur, mono):
    out = {}
    for mp, cp in cur.items():
        cur_add(out, residue_mono(mp, mono), cp)
    return out


def parity(mono):
    return sum(l[0] for l in mono) % 2


# --- V1'a: T generates translations, exactly, across parities -------
probe = [(PSI,), (DX,), (DPSI,), (D2X,), (DX, PSI), (PSI, DPSI),
         (DX, DX, PSI), ((1, 0, 0), (1, 0, 1)), ((0, 1, 1), PSI),
         ((0, 1, 0), (0, 1, 1), (1, 0, 1)), ((1, 0, 1), (1, 1, 1))]
ok = all(res_cur_mono(T, M) == {m: F(c) for m, c in d_mono(M).items()}
         for M in probe if all(l[2] == 0 for l in M))
# multi-field probes need the multi-field T
for nf in (2,):
    Tn = {}
    for j in range(nf):
        cur_add(Tn, {((0, 1, j), (0, 1, j)): F(1, 2),
                     ((1, 0, j), (1, 1, j)): F(-1, 2)})
    ok = ok and all(res_cur_mono(Tn, M)
                    == {m: F(c) for m, c in d_mono(M).items()}
                    for M in probe)
check("V1'a translation: Res(T, M) = dM exactly (mixed parities, 2 fields)", ok)

check("V1'a' Res(T, T) = dT and Res(T, G) = dG",
      residue_cur(T, T) == {m: F(c) for m, c in d_cur(T).items()}
      and residue_cur(T, G) == {m: F(c) for m, c in d_cur(G).items()})

# --- V1'b: the canary — Res(G, G) = 2T exactly ----------------------
check("V1'b canary: Res(G, G) = 2T exactly ({oint G, oint G} = 2 oint T)",
      residue_cur(G, G) == {m: 2 * c for m, c in T.items()})

# --- V1'c: SUSY variations G·psi = dX, G·dX = dpsi ------------------
check("V1'c SUSY variation: Res(G, psi) = dX, Res(G, dX) = dpsi",
      res_cur_mono(G, (PSI,)) == {(DX,): F(1)}
      and res_cur_mono(G, (DX,)) == {(DPSI,): F(1)})

# --- V1'd: fermion zero-mode bracket Res(psi, psi) = 1 --------------
check("V1'd zero mode: Res(psi, psi) = 1 (simple-pole branch)",
      residue_mono((PSI,), (PSI,)) == {(): F(1)})

# --- V1'e: Koszul basics --------------------------------------------
check("V1'e square zero + anticommutation sign",
      canonicalize((PSI, PSI)) == (None, 0)
      and canonicalize((DPSI, PSI)) == ((PSI, DPSI), -1)
      and d_mono((PSI, DPSI)) == {(PSI, (1, 2, 0)): 1})

# --- V1'f: graded antisymmetry mod d --------------------------------
# oint Res(P,Q) = -(-1)^{pq} oint Res(Q,P), i.e.
# Res(P,Q) + (-1)^{pq} Res(Q,P) is a total derivative.
pairs = [({(DX, PSI): F(1)}, {(D2X, PSI): F(1)}),          # odd, odd
         (T, {(PSI, DPSI): F(1)}),                          # even, even
         (T, G),                                            # even, odd
         ({((0, 1, 0), (1, 0, 1)): F(1)},
          {((0, 1, 1), (1, 0, 0)): F(1)})]                  # odd, odd, 2f
ok = True
for P, Q in pairs:
    pq = parity(next(iter(P))) * parity(next(iter(Q)))
    S = residue_cur(P, Q)
    cur_add(S, residue_cur(Q, P), F((-1) ** pq))
    if S:
        w2 = mono_weight2(next(iter(S)))
        nf = 1 + max(l[2] for m in S for l in m)
        ok = ok and SuperSector(w2, nf).is_total_derivative(S)
check("V1'f graded antisymmetry: Res(P,Q) + (-1)^pq Res(Q,P) in Im d", ok)

# --- V1'g: bosonic cross-check against ope_engine -------------------
def to_bos(mono):
    return tuple((m, j) for _, m, j in mono)


bos_pairs = [((DX, DX), (DX, DX)), ((DX, DX), ((0, 2, 1), (0, 1, 1))),
             (((0, 1, 0), (0, 1, 1)), ((0, 1, 0), (0, 1, 1))),
             (((0, 3, 1),), ((0, 1, 1), (0, 1, 1), (0, 2, 0))),
             ((DX, DX, (0, 2, 0)), (DX, (0, 1, 1)))]
ok = True
for P, Q in bos_pairs:
    sup = {to_bos(m): c for m, c in residue_mono(P, Q).items()}
    bos = {m: F(c) for m, c in
           ope_engine.residue_mono(to_bos(P), to_bos(Q)).items()}
    ok = ok and sup == bos
check("V1'g bosonic subalgebra matches ope_engine exactly", ok)

# --- V1'h: negative control -----------------------------------------
T_bad = {(DX, DX): F(1, 2), (PSI, DPSI): F(-2, 5)}
check("V1'h negative control: perturbed T fails translation on psi",
      res_cur_mono(T_bad, (PSI,)) != {(DPSI,): F(1)}
      and res_cur_mono(T_bad, (DX,)) == {(D2X,): F(1)})

# --- V1'i: basis generator sanity -----------------------------------
# weight2=2: dX, psi psi' pairs...  for 1 field: {(dX,), (psi dpsi)? no
# — psi dpsi has weight2 = 1+3 = 4}.  Count small cases explicitly.
b2 = gen_basis_super(2, 1)
b3 = gen_basis_super(3, 1)
b4 = gen_basis_super(4, 1)
check("V1'i basis counts (1 field, w2=2,3,4)",
      b2 == [((0, 1, 0),)]
      and set(b3) == {((1, 1, 0),), ((0, 1, 0), (1, 0, 0))}
      and set(b4) == {((0, 2, 0),), ((0, 1, 0), (0, 1, 0)),
                      ((1, 0, 0), (1, 1, 0))})

# --- V1'j: point-split normal product -------------------------------
# :TT: - 1/4 :G dG: must equal the standard quantum super-KdV I3
# density (Chistyakova-Litvinov-Orlov 2110.05870; independently
# rederived by an external audit and by this engine, 2026-08-04).
I3ps = normal_prod(T, T)
cur_add(I3ps, normal_prod(G, d_cur(G)), F(-1, 4))
I3_expected = {(DX, DX, DX, DX): F(1, 4),
               (DX, DX, PSI, DPSI): F(-3, 4),
               (DX, (0, 3, 0)): F(3, 8),
               ((0, 2, 0), (0, 2, 0)): F(-1, 4),
               (PSI, (1, 3, 0)): F(-1, 8),
               ((1, 1, 0), (1, 2, 0)): F(1, 2)}
check("V1'j point-split :TT: - 1/4 :G dG: = quantum super-KdV I3 (CLO)",
      I3ps == I3_expected)

# the correction relative to naive concatenation must NOT be d-exact
# (this is exactly the distinction the audit caught)
naive = {(DX, DX, DX, DX): F(1, 4), (DX, DX, PSI, DPSI): F(-3, 4)}
corr = dict(I3ps)
cur_add(corr, naive, F(-1))
check("V1'j' normal-ordering correction terms are not d-exact",
      bool(corr) and not SuperSector(8, 1).is_total_derivative(corr))

# graded symmetry of point-split products holds only mod d
tg = normal_prod(T, G)
cur_add(tg, normal_prod(G, T), F(-1))
check("V1'j'' :TG: - :GT: is a total derivative",
      SuperSector(7, 1).is_total_derivative(tg))

# --- V1'k: the naive-pencil ray is a genuinely different tower ------
# [I3_point-split, I3_naive] is NOT a total derivative: the legacy
# concatenation-pencil family and standard super-KdV are distinct
# commuting towers (audit finding 1; results/super_kdv_calibration.md
# section 4).
naive_I3 = {(DX, DX, DX, DX): F(1), (DX, DX, PSI, DPSI): F(-3)}
S = residue_cur(I3ps, naive_I3)
check("V1'k [I3_point-split, I3_naive] does NOT vanish mod d",
      bool(S) and not SuperSector(
          mono_weight2(next(iter(S))), 1).is_total_derivative(S))

print()
if failures:
    print(f"{len(failures)} FAILURES: {failures}")
    sys.exit(1)
print("ALL SUPER-ENGINE CHECKS PASS")
