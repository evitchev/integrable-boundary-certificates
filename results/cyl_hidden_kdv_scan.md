# The cylindrical hidden-KdV scan (post-publication item 1)

Certificate `code/cyl_hidden_kdv_scan.py` (registered; ~4 s; sweep
37/49 green, 2026-08-21).  Question, from the submitted main paper's
own caveat: could the three Z2 x O(N) cylindrical families be a
background-charge KdV / sinh-Gordon tower on a distinguished boson in
disguise -- the failure mode that fooled the companion's first version
at c = 1?  Run along the whole curve s^2 = (N-25)(N-1), both sheets.

## Setup

O(N)-invariance leaves one distinguished boson, X, hence one invariant
Feigin--Fuchs family T_Q = 1/2(dX)^2 + Q d^2X + T_Y and one invariant
sinh-Gordon pair {e^{+-aX}}.  Everything runs in the mixed
representation (X letters + O(N)-invariant Y pairs), whose basis and
total-derivative structure are N-independent, so the scan is exact at
arbitrary rational curve points -- no integer-N restriction.

## Part A -- the KdV ray: excluded at every N != 0

In-engine (component normal product at N = 1, 2; invariant projection;
reduction in the full-parity mixed sector; the two N agree):

    4 :T_Q^2: mod d = (1111)_X + (4Q^2-2)(22)_X + 2 (11)_X(11)_Y
                      + (11)(11)_Y - 2 (22)_Y   +  4Q (2)_X(11)_Y.

The last term is the point: with Y coupled, the background charge
gives the KdV ray a NONZERO X-odd class (Q d^2X dY.dY is not a total
derivative).  A Z2-even family can therefore meet the ray only at
Q = 0, and there the even-sector cross term (coefficient 2 in KdV; 6,
(11+N+s)/8, 7-N+s in Solutions 1-3) forces N = 0 on the curve.
Solver: no (lam, Q) at any sampled point with N != 0; at the formal
N = 0 point Solutions 2 and 3 ARE the c=1 KdV ray (Q = 0; positive
control), Solution 1 is not (cross term 6).

## Part B -- sinh-Gordon on X: generic exclusion, exact finite locus

First-order conservation of I_3 under e^{aX} (X contractions, Y pairs
as Taylor-shifted spectators, twisted derivative d + a dX).  Symbolic
in (N, s), each solution yields exactly TWO condition polynomials:

    Sol 1:  -6 a (a^2 - 2)                 and a * (quartic in a; N, s)
    Sol 2:  -a (a^2 - 2)(N + s + 11)       and a * (quartic)
    Sol 3:  -a (a^2 - 2)(7 - N + s)        and a * (quartic)

-- the first from the cross term (screened at a^2 = 2, or absent when
its coefficient vanishes), the second the X-part's KdV-type quartic.
Their common nonzero zeros on the curve form a FINITE set (computed
exactly; completeness is fail-closed: the ideal <c1/a, c2/a, curve> is
zero-dimensional and its lex Groebner basis is triangular, so every
point is read off and the set is asserted against the baked locus):

| Solution | points (N, s)                         | a^2     |
|----------|---------------------------------------|---------|
| 1        | (1, 0), (28, -9)                      | 2       |
| 2        | (0, 5), (28, 9), (-25/2, -45/2)       | 2       |
| 3        | (-2, -9)                              | 1 and 4 |
| 3        | (-25/2, -45/2), (0, -5), (1, 0)       | 2       |

At every generic sampled point: no nonzero root.  So the cylindrical
families are NOT hidden X-sinh-Gordon towers along the family; they
touch such structures only on a jumping locus.

Fact recorded by a control: the cross term (11)_X(11)_Y ALONE is
screened at a^2 = 2 (the double contraction; the naive single-
contraction expectation that it breaks screening is wrong).  Negative
control: a KdV ray at Q^2 = 1 plus the cross term has no common root.

## Part C -- what happens on the locus (tower level, spin 5)

* a^2 = 2 points with I_5 multiplicity one -- Sol 1 at (1,0),
  (28,-9); Sol 2 at (0,5), (28,9); Sol 3 at (0,-5), (1,0): the unique
  I_5 is screened by e^{+-sqrt2 X} as well.  Screening by the self-dual
  pair +-sqrt2 means X enters only through its stress tensor (the
  joint +-sqrt2 kernel is L(1,0)_X, the companion's Section 6): at
  these points the tower through spin 5 lives in L(1,0)_X (x)
  (Y-invariants) -- e.g. at (28, 9) the density is
  4 (7 T_X + T_Y)(T_X + T_Y) up to the usual quantum terms.  N = 1 is
  the branch point (Sol 1 = Sol 3 there, the paperclip n = -1
  reduction) and N = 0 has no Y; N = 28 is a PHYSICAL integer, with
  Solutions 1 and 2 embedded on opposite sheets.
* Sol 3 at (-2, -9) DECOUPLES: the cross term vanishes, the X-part is
  J^4 - J'^2 = the c = -2 KdV ray (kappa = -1, Q^2 = 1/4) with its
  c_{1,2} quartet {+-1, +-2} (FHW's S_2 series), the I_5 kernel has
  dimension 2 = {pure-Y I_5} + {X-KdV I_5}, fully screened at a = 1 and
  a = 2.  The family splits as KdV(c=-2) on X (+) KdV on Y.
* (-25/2, -45/2), Solutions 2 and 3: I_3 is screened at a^2 = 2, the
  I_5 kernel JUMPS to dimension 2, and -- CORRECTED 2026-08-21 after a
  Codex audit -- the e^{sqrt2 X}-screened subspace has dimension 1:
  the combination P_1 - (48/17) P_0 of the stored kernel basis
  continues the self-dual hierarchy while the other class is new.
  So this is a tower-level embedding point sitting inside an enlarged
  kernel -- the c_{1,p}-type enlargement phenomenon of the companion
  (where the KdV line persists inside a bigger kernel), here on the
  cylindrical curve.  The first version of this note reported
  screened dimension 0 ("an I_3-only coincidence"): the test had
  rescaled each current's quotient coordinates independently (numerator
  extraction per entry), which does not preserve linear dependence; the
  certificate now uses common quotient coordinates (one left-nullspace
  basis of the twisted-derivative matrix) and asserts (2, 1).
  (N = -25/2 is numerically c_{1,4}; whether that is meaningful is
  open -- recorded as a coincidence, not a claim.)
* (-2, -9) is a decoupling point for Solution 2 as well, of a different
  type: its cross term AND (dX)^4 coefficient vanish, leaving
  -2 (d2X)^2 + KdV on Y -- the free quadratic hierarchy on X
  (I_5 kernel dim 2 = pure-Y (+) pure-X).  Found by the Groebner
  completeness check as a double root at a = 0 (no screening, but a
  structural point).  Solution 1 never decouples (cross term 6).

## Verdict and consequences

1. The caveat in the submitted main paper is DISCHARGED along the
   family: no invariant KdV ray, no X-sinh-Gordon tower generically.
   The cylindrical families keep their standing as unidentified
   commuting structures -- and the paper's caution was exactly right
   to demand this check.
2. The scan produced the first exactly computed JUMPING LOCUS for a
   family of this project: finitely many curve points where the
   commuting structure degenerates into known ones (self-dual L(1,0)_X
   embeddings, a c = -2 decoupling, a multiplicity jump).  This is
   concrete input for the global jumping-loci paper (lane 7).
3. Open next steps: (a) whether the a^2 = 2 embeddings persist for
   I_7, I_9, I_11 at N = 28 (the promoted Q(t)-symbolic densities
   specialize to t = +-1/3; weight-8..12 screening in the mixed
   representation); (b) the structure of the dimension-2 kernel at
   (-25/2, -45/2); (c) the NON-invariant screening fit (O(N)-orbit
   momenta (alpha_X, beta u.Y)) -- the "towers in search of a brane"
   question proper, which needs a split-Y engine or integer-N
   component runs (N = 26, 33 are rational points).

## Audit record

Codex (2026-08-21) on commit ca18890: (1) HIGH, confirmed and fixed --
the screened-subspace test rescaled currents independently; correct
common-coordinate result at (-25/2,-45/2) is (kernel 2, screened 1)
with the combination P_1 - (48/17) P_0, reproduced in-engine before
adoption; the zero had been printed, not asserted -- now asserted.
(2) MEDIUM, confirmed and fixed in code -- public_lint's privacy screen
still exempted by basename (a manifested evil/EXPORT_RECORD.json passed
on the live v1.1.1 tag); now root-relative, regression-tested against
the reader clone (the evil file is caught); ships in v1.1.2 (verifier
defect; the exporter's build-time screen covers every staged file, so
no leak).  (3) LOW, adopted -- locus completeness via Groebner basis; a follow-up residual (the check still enumerated the triangular system with sp.solve) closed the same day: the certificate now back-substitutes explicitly -- roots read off linear factors over Q (fail-closed if any factor is nonlinear), N-elements asserted linear, a-elements asserted even in a and solved in b = a^2 -- with no root-finding solver called.

## Follow-up (2026-08-22): the N = 28 embedding persists through I_11 -- and the two solutions merge

Certificate `code/cyl_n28_tower_embedding.py` (registered; ~30 s;
Codex-audited the same day -- the audit record below).  At
(28, -9) [Solution 1] and (28, 9) [Solution 2]: I_5 (recomputed there,
unique) and the promoted Q(t)-symbolic I_7, I_9, I_11 (coordinates
A(N) + s B(N), pole-free at N = 28), each asserted to be a NONZERO
CHARGE CLASS mod d and re-verified to commute with I_3 in-engine, are
ALL conserved under e^{+-sqrt2 X}, exactly over Q(sqrt2).  Controls:
the N = 28 I_7 fails at the non-root a = sqrt3; at a generic curve
point neither I_3 nor I_7 is screened.

MERGER.  The two solutions' I_3 coincide EXACTLY at these deck-
conjugate points, and their normalized I_5, I_7, I_9, I_11 classes
coincide too: Solutions 1 and 2 are ONE family there through the
entire known range.  The common spin-3 density is exactly

    7 (dX)^4 + 6 (dX)^2 (dY.dY) + (dY.dY)^2 - 14 (d2X)^2 - 2 (d2Y.d2Y)
    = 4 [ 7 :T_X^2: + 6 T_X T_Y + :T_Y^2: ]

with no further terms -- a T-composite through and through.  (A
first version of this note wrote 4(7T_X + T_Y)(T_X + T_Y); that has
cross coefficient 8, not 6, and was wrong; the quadratic form
7u^2 + 6uv + v^2 does not factor over Q.)

Status of the statement: a FINITE-RANGE tower merger/embedding through
I_11 -- no all-spin claim.  Screening by the self-dual pair means X
enters only through its stress tensor, so through the known range the
merged N = 28 family is a commuting structure inside
L(1,0)_X (x) (O(28)-invariants of Y).  Open (recorded, not claimed):
what that structure is (not plain KdV of T_X + T_Y -- the cross-term
ratio excludes it); whether the merger and the screening persist at
all spins; whether c_Y = 28 has an independent meaning.  First work
items of the jumping-loci paper.

## Audit record (N = 28 follow-up, 2026-08-22)

Codex on commits 7571a40/3b96f37: (1) HIGH, fixed -- the certificate
tested a nonzero density, not a nonzero class (a total derivative
would pass trivially); now asserts reduce_mod_d != 0 for every
specialization (all six are genuine classes, as Codex had checked
independently).  (2) MEDIUM, fixed -- the advertised factorization was
algebraically wrong (see above).  (3) MEDIUM, adopted -- the two
N = 28 embeddings coincide; the merger is now asserted spin by spin.
(4) MEDIUM, fixed in the v2 bank -- an all-spin overstatement about
the odd Delta-classes.  (5) LOW, fixed -- I_5 now tested in this
certificate so its marker is self-contained.  (6) LOW, fixed --
explicit fail-closed flag check (survives python -O) and a
__file__-relative data path.  Wording "genuine tower-level
degeneration point" replaced by "finite-range merger/embedding".
