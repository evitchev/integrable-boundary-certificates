"""The odd-sector identification certified: [Q_{2sqrt2}, Q_{-2sqrt2}] = (sqrt2/45) J.

Codex-audit round (2026-08-20) closure: the residue identification
of the odd-sector charge J with the singlet-algebra commutator --
found by the round-2 review, previously attributed-not-certified
in odd_sector_j.py -- is here certified in-engine, together with
the exact constant (which also settles the round-3 report's
inverted normalization: J = (45/sqrt2) [Q_+, Q_-], NOT
(sqrt2/45) [Q_+, Q_-]).

The computation: with Q_{+-2sqrt2} = oint e^{+-2sqrt2 phi}, the
commutator is oint of the residue of e^{2sqrt2 phi(z)}
e^{-2sqrt2 phi(w)}.  The OPE prefactor is (z-w)^{-8}; the
normal-ordered factor is exp(2sqrt2 sum_{k>=1} (z-w)^k d^k
phi(w)/k!), so the residue current is the coefficient of
(z-w)^7 -- the complete Bell (Schur) polynomial S_7 in
t_k = 2sqrt2 d^k phi/k!.  Checks (fail-closed):

 1. The parity-EVEN part of the residue (even part counts of the
    partitions of 7) is a total derivative -- forced by parity,
    since the commutator charge is Z2-odd.
 2. The parity-ODD part is nonzero mod d and EXACTLY proportional
    to J, with stripped ratio 1/45:
       [Q_{2sqrt2}, Q_{-2sqrt2}] == (sqrt2/45) J  mod d.
 3. Negative control: the perturbed J' of odd_sector_j.py
    (64/7 -> 65/7) is NOT proportional to the residue class.

Together with odd_sector_j.py (J nonzero, commutes with the tower,
conservation pattern) and the Jacobi identity, this makes the
all-spin statement self-contained: J commutes with every charge
annihilated by both +-2sqrt2 conditions.
"""
import sys
import math
from collections import Counter
from fractions import Fraction as F

from families import L
from ope_engine import Sector

failures = []


def check(name, ok):
    print(f"{'PASS' if ok else 'FAIL'}  {name}", flush=True)
    if not ok:
        failures.append(name)


def partitions(n, mx=None):
    if mx is None:
        mx = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, mx), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest


# residue current split by parity of the part count j:
# (2sqrt2)^j = 2^{3j/2} (j even) or sqrt2 * 2^{(3j-1)/2} (j odd)
C_even, C_odd = {}, {}
for lam in partitions(7):
    j = len(lam)
    m = Counter(lam)
    denom = 1
    for k, mk in m.items():
        denom *= (math.factorial(k) ** mk) * math.factorial(mk)
    mono = L(*[(k, 0) for k in lam])
    tgt = C_even if j % 2 == 0 else C_odd
    num = 2 ** (3 * j // 2) if j % 2 == 0 else 2 ** ((3 * j - 1) // 2)
    tgt[mono] = tgt.get(mono, F(0)) + F(num, denom)

J = {L(*((1, 0),) * 7): F(64, 7),
     L((1, 0), (1, 0), (1, 0), (2, 0), (2, 0)): F(-40),
     L((1, 0), (3, 0), (3, 0)): F(1)}

sec = Sector(7, 1, False)
ve = sec.reduce_mod_d(sec.vec(C_even))
vo = sec.reduce_mod_d(sec.vec(C_odd))
vJ = sec.reduce_mod_d(sec.vec(J))

# ---- 1. parity kills the even part ------------------------------------
check("parity-even part of the residue is a total derivative",
      not any(ve))

# ---- 2. odd part == (1/45) J (sqrt2 stripped) -------------------------
check("parity-odd part nonzero mod d", any(vo))
check("J nonzero mod d", any(vJ))
iJ = next(i for i, x in enumerate(vJ) if x)
ok = vo[iJ] != 0
if ok:
    r = F(vo[iJ], vJ[iJ])
    ok = (r == F(1, 45)
          and all(vo[i] == r * vJ[i] for i in range(len(vJ))))
check("[Q_{2sqrt2}, Q_{-2sqrt2}] == (sqrt2/45) J mod d "
      "(stripped ratio exactly 1/45)", ok)

# ---- 3. negative control ----------------------------------------------
Jbad = dict(J)
Jbad[L(*((1, 0),) * 7)] = F(65, 7)
vb = sec.reduce_mod_d(sec.vec(Jbad))
okneg = vb[iJ] != 0
if okneg:
    rb = F(vo[iJ], vb[iJ])
    okneg = not all(vo[i] == rb * vb[i] for i in range(len(vb)))
check("negative control: perturbed J' is not proportional to the "
      "residue class", okneg)

print()
if failures:
    print(f"J-COMMUTATOR IDENTIFICATION CHECKS FAILED ({len(failures)}):")
    for f in failures:
        print(" -", f)
    sys.exit(1)
print("J-COMMUTATOR IDENTIFICATION CERTIFIED "
      "([Q_{2sqrt2}, Q_{-2sqrt2}] = (sqrt2/45) J mod d; "
      "even part total derivative; constant exact)")
