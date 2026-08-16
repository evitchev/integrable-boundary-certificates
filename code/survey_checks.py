"""Computational corroborations for the survey section of the draft.

A. Commuting-variety structure: for one free boson (Z2 sector), the
   variety of commuting pairs of spin-3 charges in Q^3 x Q^3 is exactly
   the rank-one (proportionality) determinantal locus: the commutator
   class equals (ab' - a'b) K with K a fixed nonzero spin-6 charge
   class.
B. Arc-space/Hilbert-series check: machine dims of B^s and Q^sigma for
   the Z2 sector of one boson match the closed-form character
   chi_e(q) = ( prod 1/(1-q^k) + prod 1/(1+q^k) ) / 2  (even number of
   letters), with dim Q^sigma = [q^(sigma+1)]chi_e - [q^sigma]chi_e.
C. Resurgence of the Wronskian constants: Vitchev's beta_{2k-1} =
   B_{2k}/(2k(2k-1)2^(2k-1)) (his eq. (72)) satisfy
   (16 pi^2) |beta_{2k+1}/beta_{2k-1}| / ((2k)(2k-1)) -> 1,
   i.e. Gevrey-1 growth with nearest Borel singularity at |t| = 4 pi
   in the 1/kappa plane -- the Stokes structure of log Gamma(2 kappa).
"""

from fractions import Fraction as F

from ope_engine import Sector, gen_basis, residue_cur

# ---------------------------------------------------------------- A
print("== A. spin-(3,3) commuting variety, one boson, Z2 sector ==")
R1 = {((2, 0), (2, 0)): F(1)}
R2 = {((1, 0), (1, 0), (1, 0), (1, 0)): F(1)}
sec4 = Sector(4, 1, z2=True)
assert sec4.q_dim() == 2
sec7 = Sector(7, 1, z2=True)
K11 = residue_cur(R1, R1)
K22 = residue_cur(R2, R2)
K12 = residue_cur(R1, R2)
K21 = residue_cur(R2, R1)
print("  [R1,R1], [R2,R2] total derivatives:",
      sec7.is_total_derivative(K11), sec7.is_total_derivative(K22))
S = dict(K12)
for m, c in K21.items():
    S[m] = S.get(m, 0) + c
print("  K12 + K21 total derivative (antisymmetry):",
      sec7.is_total_derivative({m: c for m, c in S.items() if c}))
print("  K12 itself a total derivative:", sec7.is_total_derivative(K12))
print("  => [aR1+bR2, a'R1+b'R2] = (ab'-a'b) K, K != 0:")
print("     commuting variety = {ab' = a'b}, the rank-1 determinantal locus")

# ---------------------------------------------------------------- B
print("\n== B. character / Hilbert-series check, one boson, Z2 sector ==")
TOP = 17
# chi_e coefficients via exact polynomial products
prod_m = [F(0)] * (TOP + 1)
prod_m[0] = F(1)
prod_p = list(prod_m)
for k in range(1, TOP + 1):
    # multiply prod_m by 1/(1-q^k), prod_p by 1/(1+q^k), truncating
    for s in range(k, TOP + 1):
        prod_m[s] += prod_m[s - k]
    new_p = list(prod_p)
    for s in range(k, TOP + 1):
        new_p[s] = prod_p[s] - new_p[s - k]
    prod_p = new_p
chi = [(a + b) / 2 for a, b in zip(prod_m, prod_p)]
okB = okQ = True
for s in range(2, TOP):
    if len(gen_basis(s, 1, z2=True)) != chi[s]:
        okB = False
for sigma in range(2, 13):
    if Sector(sigma + 1, 1, z2=True).q_dim() != chi[sigma + 1] - chi[sigma]:
        okQ = False
print(f"  dim B^s == [q^s] chi_e for s = 2..{TOP-1}:", okB)
print("  dim Q^sigma == [q^(sigma+1)]chi_e - [q^sigma]chi_e "
      "for sigma = 2..12:", okQ)

# ---------------------------------------------------------------- C
print("\n== C. Borel growth of the Wronskian constants beta_{2k-1} ==")
import sympy as sp
pi2_16 = 16 * sp.pi**2


def beta(k):
    return sp.bernoulli(2 * k) / (2 * k * (2 * k - 1) * 2 ** (2 * k - 1))


for k in (5, 10, 20, 40, 80):
    ratio = abs(beta(k + 1) / beta(k))
    val = sp.N(pi2_16 * ratio / ((2 * k) * (2 * k - 1)), 12)
    print(f"  k = {k:3d}:  (16 pi^2)|b_(2k+1)/b_(2k-1)| / (2k)(2k-1) = {val}")
print("  -> limit 1: Gevrey-1, nearest Borel singularity |t| = 4 pi")
