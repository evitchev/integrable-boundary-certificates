# NZERO1 Prediction (sealed before computation)

Date: 2026-10-07, Qwen seat.

## Mechanism (from engine trace)

At t=4 (k=-6, n=-9) spin 9 and at t=7 (k=-3, n=-3) spins 3, 9, the WKB engine's
integrate_reg() step produces a raw charge containing ONLY PX^0 (PI^2-only) monomials.
The monic normalisation in score() divides by the PX^(spin-1) coefficient, which is zero.
This is NOT the Beta prefactor being infinity/infinity in a way that kills the charge;
rather, the pole-residue replacement in integrate_reg combined with the specific Gamma(-c)
values at these parameters causes all top-power terms to have zero coefficient.

The file:line of the zero:
- thwkb.py line ~218-225 (integrate_reg): the pole-residue replacement
- ext_fibres.py line ~47 (score): gnorm(ch[s_]) == 0 check -> "WKB charge identically 0"
- The actual vanishing: the sum over poly items in run() (ext_fibres.py line ~40-42)
  produces a polynomial with no PX^(spin-1) term.

## Prediction for the exact limit

The brief asks: compute the coefficient as an exact LIMIT along the family t->4 (spin 9)
and t->7 (spins 3, 9), and compare with the certified charge.

Prediction: The limit ALONG THE FAMILY (t = t0 + eps, eps -> 0) of the WKB charge,
properly normalised, will RECOVER the certified charge. The reason: at t = t0 +/- 1/100
the agreement is exact (per the brief's observation), so the function is continuous through
the exceptional point. The "zero" is an artifact of the discrete monic normalisation at the
exact rational point where the Gamma functions have coincident poles; approaching it
generically gives the correct finite value.

Specifically:
- At t=4, spin 9: limit_{t->4} [monic(WKB_charge(t))] = certified(2, 4)[9]
- At t=7, spin 3: limit_{t->7} [monic(WKB_charge(t))] = certified(2, 7)[3]
- At t=7, spin 9: limit_{t->7} [monic(WKB_charge(t))] = certified(2, 7)[9]

The mechanism: the Beta prefactor B(-k*nu, nu*(1+k)) at the exceptional point has both
Gamma factors diverging, but their RATIO (which is what enters the charge after monic
normalisation) has a finite limit. This is the standard "pole/pole = finite" mechanism
from the Gamma function: Gamma(x+eps)/Gamma(x) -> 1/Res(Gamma, x) as eps -> 0... 
actually more precisely, the product of residues from the coincident poles cancels
the leading divergence, leaving a finite value equal to the directional derivative
of the regularised expression.

## Controls (predicted outcomes)

1. Regular cell (t=4, spin 7): direct WKB value should match certified(2,4)[7] EXACTLY
   (this is already in the run log; confirming as control).
2. Non-exceptional neighbour (t=4+1/100): direct WKB value matches certified(2, 4.01)[...]
   at all spins (already observed per brief).
3. Tamper (wrong k in Beta prefactor): if we replace k with k+1/10 in the prefactor,
   the limit should BREAK (not equal certified).

## What this means for the paper

The "Exceptional cells" paragraph in Sec. 6 of the paper says these cells have
"nu = spin/n a negative integer and k*nu a positive integer; the common Beta prefactor
B(-k nu, nu(1+k)) is of the form infinity/infinity." The correct statement is:

The WKB coefficient at these points is the FINITE LIMIT of the generic-t expression.
The engine's monic normalisation fails at exactly these rational points because it
divides by a vanishing leading coefficient, but the charge itself has a well-defined
finite value obtainable as the t -> t0 limit. This is analogous to the nZ (multiples of n)
cells where the operator vanishes and the charge is obtained as a directional limit.

