# Hand-off spec: independent Mathematica replication (verification triangle)

For the clean-context implementer (cognitive's item-4 harness).  This
document (shipped as wolfram/HANDOFF-SPEC.md), the two submitted
papers, hep-th/0404195, and the target list `code/wolfram_targets.json`
are the ONLY inputs.

## Independence rules (non-negotiable)

1. Do NOT read anything under `code/` of integrable-boundary, nor
   `results/*.json`, nor the independent/ engine.  The point of this
   exercise is an unrelated implementation, vendor, and algorithm
   lineage.  If a convention below is ambiguous, REPORT the ambiguity
   in the output (field `notes`) -- do not resolve it by looking at
   our code.
2. Preferred engine: Thielemans' OPEdefs (check its license and that
   it loads under the installed Mathematica / Wolfram Engine first).
   Fallback: a from-scratch engine written from the conventions here.
3. Execute at every step.  A script that has never run is not a
   result.  Record engine version, package version, script SHA-256,
   and a run timestamp in the output.
4. Output: a single JSON file `results/wolfram_triangle.json` with one
   entry per claim id (schema below).  Our runner
   (`code/wolfram_triangle.py`) compares it fail-closed against the
   baked targets.

## Conventions (from the papers)

* One free boson phi with <d phi(z) d phi(w)> = (z-w)^{-2} (so
  <d^m phi d^n phi> = (-1)^{m-1} (m+n-1)! (z-w)^{-(m+n)}).  Several
  bosons X, Y_1..Y_N: independent, each with this normalization.
* A CHARGE is oint P for a differential polynomial P; charges are
  identified modulo total derivatives: the CLASS of P is P mod Im d.
  "A = B mod d" means A - B is a total derivative.  Spin sigma charge
  = density of weight sigma + 1.
* Normal product (point-split, the OPEdefs NO[A,B] convention):
  :AB:(w) = oint dz/(z-w) A(z) B(w), i.e. the (z-w)^0 term of the OPE
  A(z)B(w).  Nested products are right-to-left as written:
  :T:TT:: = NO[T, NO[T, T]].
* Feigin--Fuchs stress tensor: T_Q = 1/2 (d phi)^2 + Q d^2 phi,
  central charge c(Q) = 1 - 12 Q^2.  Vertex operator V_alpha = e^{alpha
  phi}, weight h = alpha(alpha - 2Q)/2; at Q = 0, h = alpha^2/2.
* Commutator of charges: [oint A, oint B] = oint_w Res_{z->w} A(z)B(w)
  (the single-pole term); two charges commute iff that residue
  density is a total derivative.
* First-order conservation under a perturbation V_alpha (Zamolodchikov's
  criterion, as used in the companion Section 4): oint P is conserved
  iff Res_{z->w} P(z) V_alpha(w) = (d + alpha d phi) J (w) V_alpha(w)
  for some differential polynomial J -- a TWISTED total derivative.
  (Multi-boson: alpha.dX in the twist.)
* Cylindrical families (hep-th/0404195 Sec. 3.3, eqs. (53), (55),
  (57), with the normalization rule stated in the main paper's
  Section "The cylindrical families at general N": the Y-part reads
  (dY.dY)^2 - 2 (d2Y.d2Y), i.e. "4 T_Y^2" is to be understood as
  that expression, and the X-dependent terms carry the coefficients
  below).  Fields X, Y_1..Y_N; dY.dY = sum_i (dY_i)^2 etc.  On the curve
  s^2 = (N - 25)(N - 1):
    Sol 1: P = (2+N+s)/3 (dX)^4 + 6 (dX)^2 (dY.dY) + (dY.dY)^2
               + (-9-3N-s)/6 (d2X)^2 - 2 (d2Y.d2Y)
    Sol 2: P = (117+2N+N^2+(15+N)s)/192 (dX)^4 + (11+N+s)/8 (dX)^2 (dY.dY)
               + (dY.dY)^2 - (409+46N+N^2+(N-5)s)/192 (d2X)^2 - 2 (d2Y.d2Y)
    Sol 3: P = (dX)^4 + (7-N+s) (dX)^2 (dY.dY) + (dY.dY)^2
               + (13-27N+2N^2+(5-2N)s)/6 (d2X)^2 - 2 (d2Y.d2Y)
  The targets below only need these at (N, s) = (28, -9) and (28, 9).
  Permitted simplification: every target involving Y is a statement
  about O(N)-invariant monomials whose verification does not depend on
  the number of explicit Y bosons (no Y contraction with e^{aX}; the
  :T_Y^2: quantum correction is N-independent), so implement with
  N_explicit = 2 and, as a check, N_explicit = 3, keeping the N = 28
  coefficient VALUES.

## Targets (first tier; each a boolean, a number, or a short list)

T1  [boolean]  :T_Q T_Q: = 1/4 [ (d phi)^4 + (4Q^2 - 2)(d2 phi)^2 ]  mod d,
    identically in Q (symbolic Q).           (companion, Theorem 1(i))
T2  [boolean]  at Q^2 = 9/8:
    :T_Q :T_Q T_Q:: + (7/8) :dT_Q dT_Q: = (1/8) P6  mod d,
    P6 = (d phi)^6 + (5/2)(d phi)^2 (d2 phi)^2 + (31/24)(d3 phi)^2.
                                               (companion, Theorem 1(ii))
T3  [exact number]  with Q = 0 and a = 2 sqrt2:
    Res_{z->w} e^{a phi(z)} e^{-a phi(w)}  =  r * J_density  mod d,
    J_density = (64/7)(d phi)^7 - 40 (d phi)^3 (d2 phi)^2 + d phi (d3 phi)^2.
    Target: r = sqrt(2)/45.                    (companion, Section 3)
T4  [set]  the first-order conservation condition of
    oint[(d phi)^4 + (5/2)(d2 phi)^2] under e^{alpha phi} (Q = 0) holds
    exactly for alpha^2 in {1/2, 8} (and alpha = 0).
                                               (companion, eq. (cond))
T5  [list of 13 integers]  dimension, at each weight w = 2..14, of the
    space of charge classes (all differential polynomials in one boson,
    any parity, mod d) conserved to first order under BOTH e^{+sqrt2
    phi} and e^{-sqrt2 phi}:  [1,0,1,0,2,0,3,1,4,2,7,3,10].
                                               (companion, Section 6)
T6a [boolean]  Sol 1 at (N,s) = (28,-9) and Sol 2 at (28, 9) give the
    SAME density (coefficient by coefficient).
T6b [boolean]  that density equals 4 [ 7 :T_X T_X: + 6 T_X T_Y + :T_Y T_Y: ]
    mod d, with T_X = 1/2 (dX)^2, T_Y = 1/2 dY.dY (point-split products;
    T_X T_Y has no contractions).
T6c [boolean pair]  that density is conserved to first order under
    e^{+-sqrt2 X} (true) and NOT under e^{+-sqrt3 X} (false).
                                      (results/cyl_hidden_kdv_scan.md)

## Output schema (per claim id)

  {"T1": {"value": true, "exact": true, "notes": "",
          "engine": "Mathematica 14.1 / OPEdefs 3.x",
          "script_sha256": "...", "run": "2026-..."}, ...}

Numbers as exact strings ("Sqrt[2]/45"); lists as JSON arrays.  Any
convention ambiguity goes in `notes`, never silently resolved.

## Conventions pinned after the first replication (2026-08-22)

The implementer recorded three choices the spec had left open; all
agree with the targets and are now the spec's definitions:
* T3: the residue is Res_{z->w} e^{+a phi(z)} e^{-a phi(w)} in that
  literal order (the opposite order flips the sign of r).
* T4: the value is the set of alpha^2 over the NONZERO roots; alpha = 0
  (the trivial vertex operator) is always a root and is excluded.
* T6c: each boolean means "conserved under BOTH signs of the
  momentum"; the density is even in X, so the two signs agree.

Outcome of the first replication: Wolfram Engine 15.0 + OPEdefs 3.1
(beta 4), all eight targets agree (WOLFRAM TRIANGLE VERIFIED; runner
re-execution 61 s).  Follow-up requested of the implementer: record
the script path in the output JSON repository-relative, not absolute
(the absolute local path is a privacy-screen hit; the exporter would
redact it at build time, but the artifact should be clean at source).
