# Q(t)-promotion of sampled existence results (2026-07-29)

All eleven target charges are now constructed as explicit symbolic
densities over the rationalized double cover
$s^2 = (N-25)(N-1)$, $N = \frac{t^2-25}{t^2-1}$, $s = \frac{-24t}{t^2-1}$,
and the defining commutators verified **identically in $t$** — proving
existence at the generic point of the parameter curve. Combined with
the sampled kernel dimensions $\le 1$ (rank semicontinuity), the
generic kernels are **exactly one-dimensional**: existence and
uniqueness both hold generically.

Pipeline (`code/promote_qt.py`, 11 parallel workers, 160 min wall /
660 min CPU): exact kernel solves at 31 pairs $\pm t$ (62 points per
charge; points with kernel dimension $\neq 1$, e.g. the cylindrical
decoupling point, are skipped automatically); canonical density
coordinates split as $A(N) + s\,B(N)$ via the deck map $t \to -t$;
$A, B$ reconstructed as rational functions of $N$ with held-out
verification; final symbolic commutator reduction modulo total
derivatives over $\mathbb{Q}(t)$.

| charge | verified | density coordinates |
|---|---|---|
| A₂⁽²⁾ I₁₁ (spin-12 density) | identically in t | 31 |
| A₂⁽²⁾ I₁₃ (spin-14 density) | identically in t | 59 |
| cyl Sol 1/2/3 I₇ (spin-8) | identically in t | 29 each |
| cyl Sol 1/2/3 I₉ (spin-10) | identically in t | 68 each |
| cyl Sol 1/2/3 I₁₁ (spin-12) | identically in t | 160 each |

Exact coefficient data ($A_i(N), B_i(N)$ per canonical coordinate,
normalization coordinate recorded): `qt_promotion.json` (~290 KB).

Consequences propagated to the draft:
- the $A_2^{(2)}$ spin content through 13 and the cylindrical all-odd
  towers through 11 are now **generic-point theorems**, not sampled
  observations;
- with the symbolic densities available, the outstanding cylindrical
  mutual-commutativity pairs ($[I_5,I_{11}]$, $[I_7,I_9]$,
  $[I_7,I_{11}]$, $[I_9,I_{11}]$) become finite symbolic checks — the
  natural next computation.
