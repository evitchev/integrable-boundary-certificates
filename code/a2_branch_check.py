"""Independent verification of the branch-point degeneration claim from a
parallel investigation: at N = 1 (c = 1, where sqrt((N-25)(N-1)) = 0) the
A_2^(2)/Phi_{1,2} family (spins 1, 6n±1) allegedly ENHANCES to a full
odd-spin tower: a spin-3 charge I_3' (distinct from KdV I_3) and a
spin-9 charge appear, filling the 6n±1 gaps.

We verify with our independently written engine:
  1. the A_2^(2) I_5 at N=1 admits genuine spin-3 charge(s) I_3';
  2. I_3' is NOT the KdV I_3 (mod total derivatives);
  3. the unique spin-7 charge I_7 commuting with I_5 satisfies
     [I_3', I_7] = 0 (the tower is coherent);
  4. a genuine spin-9 charge commuting with I_5 exists (gap 9 fills in);
  5. spectral scan of ad_{I_5^{A2}} and ad_{I_3'} across spins, compared
     against the KdV ad_{I_3} scan (two distinct all-odd towers at c=1).

A_2^(2) I_5 density, eq. (34) of hep-th/0404195 at N=1 (sq = 0):
  (1/96)(107+17N) (33) - 6 (22)(11) - (1/8)(-85+17N) (12)(12) + (11)^3
with the N=1 monomial collision (22)(11) == (12)(12).
"""

from fractions import Fraction as F

from ope_engine import Sector, gen_basis, residue_cur, rref, fmt_cur
from families import L, KDV_P4
from solve_charges import genuine_kernel


def acc(pairs):
    out = {}
    for mono, c in pairs:
        out[mono] = out.get(mono, 0) + c
    return {m: c for m, c in out.items() if c}


N = 1
A2_P6 = acc([
    (L((3, 0), (3, 0)), F(107 + 17 * N, 96)),
    (L((2, 0), (2, 0), (1, 0), (1, 0)), F(-6)),
    (L((1, 0), (2, 0), (1, 0), (2, 0)), -F(-85 + 17 * N, 8)),
    (L((1, 0),) * 6, F(1)),
])
print("A2^(2) I_5 density at N=1:")
print(fmt_cur(A2_P6, nfields=1))


def kernel_at(P_low, sigma):
    cand = gen_basis(sigma + 1, 1, z2=True)
    if not cand:
        return []
    cs = Sector(sigma + 1, 1, z2=True)
    s_low = sum(m for m, _ in next(iter(P_low)))
    tgt = Sector(s_low + sigma + 1 - 1, 1, z2=True)
    return genuine_kernel(P_low, cand, tgt, cs)


failures = []


def check(name, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {name}", flush=True)
    if not ok:
        failures.append(name)


def same_charge(c1, c2, s):
    sec = Sector(s, 1, z2=True)
    v1 = sec.reduce_mod_d(sec.vec(c1))
    v2 = sec.reduce_mod_d(sec.vec(c2))
    rows = [list(v1), list(v2)]
    _, piv = rref(rows, len(v1))
    return len(piv) < 2  # proportional mod d


# 1-2: spin-3 charges commuting with I_5^{A2}
sols3 = kernel_at(A2_P6, 3)
check("multiplicity one at sigma=3", len(sols3) == 1)
for s in sols3:
    print(fmt_cur(s, nfields=1))
    check("I_3' is NOT the KdV I_3 mod d",
          not same_charge(s, KDV_P4, 4))

# 3: I_7 and coherence
sols7 = kernel_at(A2_P6, 7)
check("multiplicity one at sigma=7", len(sols7) == 1)
if sols3 and sols7:
    I3p, I7 = sols3[0], sols7[0]
    res = residue_cur(I3p, I7)
    check("[I_3', I_7] = 0",
          Sector(4 + 8 - 1, 1, z2=True).is_total_derivative(res))

# 4: spin-9
sols9 = kernel_at(A2_P6, 9)
check("multiplicity one at sigma=9", len(sols9) == 1)

# extra coherence: [I_3', I_9], [I_7, I_9]
if sols3 and sols9:
    check("[I_3', I_9] = 0",
          Sector(4 + 10 - 1, 1, z2=True).is_total_derivative(
              residue_cur(sols3[0], sols9[0])))
if sols7 and sols9:
    check("[I_7, I_9] = 0",
          Sector(8 + 10 - 1, 1, z2=True).is_total_derivative(
              residue_cur(sols7[0], sols9[0])))

# also: does I_3' commute with KdV I_3? and cross-tower separation
if sols3:
    check("[I_3', I_3^KdV] != 0 (distinct towers)",
          not Sector(4 + 4 - 1, 1, z2=True).is_total_derivative(
              residue_cur(sols3[0], KDV_P4)))
if sols7:
    check("[I_3^KdV, I_7^A2] != 0 (distinct towers)",
          not Sector(4 + 8 - 1, 1, z2=True).is_total_derivative(
              residue_cur(KDV_P4, sols7[0])))

# 5: spectral scans
print("\n(5) spectral scans (genuine kernel dims):")
print("sigma                    | 2 3 4 5 6 7 8 9 10")
row = [len(kernel_at(A2_P6, s)) for s in range(2, 11)]
print("ker ad_{I_5^A2} at N=1   |", " ".join(str(x) for x in row))
check("all-odd multiplicity-one pattern sigma=2..10 (0,1 alternating)",
      row == [0, 1, 0, 1, 0, 1, 0, 1, 0])
if sols3:
    row = [len(kernel_at(sols3[0], s)) for s in range(2, 11)]
    print("ker ad_{I_3'}  at N=1    |", " ".join(str(x) for x in row))
row = [len(kernel_at(KDV_P4, s)) for s in range(2, 11)]
print("ker ad_{I_3^KdV} (ref)   |", " ".join(str(x) for x in row))

print()
if failures:
    print(f"FAILURES: {failures}")
    print("A2 BRANCH MECHANISM FAILED")
    import sys
    sys.exit(1)
print("A2 BRANCH MECHANISM COMPLETE (fail-closed)")
