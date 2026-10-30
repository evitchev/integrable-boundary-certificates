# q8_density_test: the q = 8 test of the singlet reading at density level (t = 9/7, weight 14)

## What the computation shows

At q = 8 (beta^2 = 1/8, t = 9/7) the extra sheet-3 class has weight 14. The computation writes the two bosons X and phi as
V = X + phi and Z = phi - X, two independent Coulomb gases, each at c = -25/2 = c_(1,4) with <V,V> = <Z,Z> = 2, and finds, in exact
rational arithmetic (Python fractions, sympy rationals, python-flint fmpq matrices; no floating point):
the kernel K_14 of ad I_3 on weight-14 fields modulo total derivatives has dimension 2 and contains the class whose vacuum
eigenvalue is the certified Solution-3 weight-14 record at t = 9/7; among the 225 products f(V) g(Z) of a weight-7 field of V with
a weight-7 field of Z, those in the kernel modulo total derivatives form a one-dimensional space, spanned by H_V H_Z, where
H = Res e^{2V}(z) e^{-2V}(w) is the weight-7 singlet generator; the vacuum eigenvalue of H_V is binom(2 P_V + 3, 7) exactly, and
that of H_V H_Z is -(64/99225) prod_{j=-3}^{3} (P_X + pi - j/2)(P_X - pi - j/2) exactly.

Notation of the paper: U, W = (phi +- X)/sqrt2, so V = sqrt2 U and Z = sqrt2 W; H_V H_Z here is H_U H_W there; P_V = P_X + pi = P_+
and P_Z = pi - P_X = -P_-, so binom(2 P_Z + 3, 7) = binom(3 - 2 P_-, 7).

## Statements of the paper that rest on it

Jumping-loci paper, Section 4 (The fusion-product law), the paragraph "At q=8 (t=9/7, weight 14) Conjecture 1 was tested against
predictions recorded before the computation": "The kernel of ad I_3 modulo total derivatives has dimension 2, as in the records,
and contains the class of the certified spin-13 record; the span of the 225 products of a weight-7 field of U with a weight-7 field
of W meets the kernel modulo total derivatives in a one-dimensional subspace, spanned by H_UH_W. The factors have the eigenvalues
binom(2P_+ + 3, 7) and binom(3 - 2P_-, 7) [...], and the eigenvalue of H_UH_W is -(64/99225) prod_{j=-3}^{3}(P_X+pi-j/2)(P_X-pi-j/2),
all exactly [lab result]. The conjecture therefore holds at the density level at q=4,6,8." The abstract's "at q=4,6,8 the extra
class is the product of two such generators" rests on it for q = 8.

Where each part is in the records:

- dimension 2: `records/w3_kernel.log`, GATE G2 (rank 893 of 895 over Q and modulo two primes); `records/w3_kernel.json`, "dim_K14".
- the certified spin-13 class in the kernel: `records/w3_kernel.log`, GATE G3.
- the 225 products: `records/w2_bilinear.log`, SEALED TEST P2 (rank 224 of 225; the relation vector equals h_i h_j, h = the coefficients of H).
- the factor eigenvalue: `records/w1_gates.log`, SEALED TEST P1.
- the eigenvalue of H_V H_Z: `records/w2_bilinear.log`, SEALED TEST P3; the identity binom(2P_V+3,7) binom(2P_Z+3,7) = -(64/99225) E_3
  also in `records/s1_frame2.log`.

## Files

- `README.md`: this file.
- `run.sh`: re-runs the five scripts into a new directory.
- `compare.sh`: compares the new outputs with `records/`.
- `code/s1_frame2.py`: step 1, the frame at t = 9/7 (V, Z, kappa = 2, screenings, stress tensors, momenta) derived from the
  screening-frame record, and the binomial identities; exact sympy algebra, no kernel.
- `code/w0_engine.py`: engine controls: the engine and the Euler-criterion route reproduce the q = 4 (kappa = 1) and q = 6
  (kappa = 3/2) results; the complement basis C_14 (895 monomials) is checked; the Euler criterion agrees with the d-echelon.
- `code/w1_gates.py`: the kappa = 2 frame by OPE, gate G1 (I_3), H and its properties, prediction P1, P4a.
- `code/w2_bilinear.py`: predictions P2 (the 225 products), P3 (H_V H_Z in K_14, its eigenvalue), P4b; controls C2-C6.
- `code/w3_kernel.py`: the full kernel K_14 on C_14, gates G2 and G3, cross-checks X1-X5; writes `w3_kernel.json`.
- `code/q8_lib.py`: shared library: the cylinder vacuum-eigenvalue rule, the Euler criterion for total derivatives, the
  complement basis, exact linear algebra with python-flint, readers of the repository records.
- `code/vendor/kwick.py`, `code/vendor/wick.py`: the free-field Wick engine (byte copies); `code/vendor/SOURCE_SHA256.txt`: their
  sha256 and source.
- `inputs/z2_kernel_t9_7_K7.json`: the two kernel classes at t = 9/7 computed in another frame (one boson X and a Virasoro algebra
  at central charge N), read for cross-check X1 only; `inputs/SOURCE_SHA256.txt`: its sha256 and source.
- `records/s1_frame2.log`, `records/w0_engine.log`, `records/w1_gates.log`, `records/w2_bilinear.log`, `records/w3_kernel.log`:
  the logs of the archived run (October 2026), each ending with the exit code ("EXIT n").
- `records/w3_kernel.json`: sizes, ranks, dim K_14, the eigenvalues of the two kernel classes, the certified weight-14 eigenvalue at
  t = 9/7, the eigenvalue of H_V H_Z, and the two kernel representatives in C_14 coordinates.

Repository files read by the scripts (their sha256 is checked by `run.sh` before anything runs):
`results/lab/first_order/vev_w4_points.json` (gate G1 and an engine control), `results/lab/vev/vev_sol3_w10.json` (an engine
control at q = 6), `results/lab/anchor13/vev_profile_sol3_w14.json` (gate G3).

## Reading the logs

- SEALED TEST marks a prediction recorded before the computation (P1, P2, P3, P4a, P4b); GATE, a condition required for the test to
  be meaningful (G1-G3); CONTROL, an engine check or a deliberately altered input that must be detected; CROSS-CHECK, a reported
  comparison that decides nothing.
- A field is printed as {monomial: coefficient}; a monomial is a sorted tuple of letters (c, q) = d^{q+1} Phi_c with c = 1 for V and
  c = 2 for Z; P_A and P_B in printed eigenvalues are P_V and P_Z.
- WUW1 and WUW2 (the q = 4 and q = 6 cases), VIR2 (the screening frame) and VIR5d (the weight-14 kernel in the other frame), and the
  paths results/lab/large_n/... in comments, name earlier archived computations of the same project; their records are in this
  repository under opers-2026-10/tree/results/lab/large_n/.
- Exit codes: every step exits 0 except `w2_bilinear.py`, which exits 1, as recorded. Its control C4 adds (dV)^2 d^2V d^3V/10 to H.
  The altered product leaves the kernel, as required, but the eigenvalue tests cannot see the change: that monomial's vacuum
  eigenvalue is identically 0 (writing it as a product of d^{q+1}V, the q >= 1 parts add up to 3, and every monomial whose q >= 1
  parts add up to an odd number has eigenvalue 0). The control was therefore vacuous for P1 and P3 and is reported as failed; every
  prediction in w2 holds. `run.sh` expects exit 1 for this step.
- The shipped logs differ from the archived ones in two places, and the scripts were changed to print the same: the line after
  line 1, which printed the sha256 of a file that is not part of this folder, is removed (w0-w3); in `w1_gates.log` one internal
  label is replaced by "predicted".

## How to re-run

From the repository root (any new output directory outside this folder):

    bash results/lab/jumping_loci_records/q8_density_test/run.sh ../q8_out
    bash results/lab/jumping_loci_records/q8_density_test/compare.sh ../q8_out

Set `PYTHON=/path/to/python3` to choose the interpreter (default `python3`; it is called as `$PYTHON -I -B`).

Environment of the archived run and of the re-run: Python 3.12.8; sympy 1.14.0, mpmath 1.3.0, python-flint 0.9.0 (FLINT 3.6.0);
numpy (2.5.2 installed) is not used. python-flint is required here (the repository's requirements.txt lists it as optional).
Also needed: bash, GNU coreutils (timeout, sha256sum), GNU sed; Linux.

Runtime and memory: about 4 minutes on one core, one step at a time (a re-run on 2026-10-30 took 3 min 58 s: s1_frame2 a few
seconds, w0 about 15 s, w1 under 1 s, w2 about 100 s, w3 about 110 s); peak resident memory about 1 GB (0.97 GiB, in the largest
step). Each step runs under ulimit -v 8000000 KiB with a timeout of 30 to 120 minutes.

## What "reproduced" means

`run.sh` ends with "every step exited as recorded" (each exit code equals the recorded one) and `compare.sh` prints
"ALL IDENTICAL". `w3_kernel.json` is compared byte for byte. The logs are compared after the same mask on both sides:

1. line 1 of each log, "<script>.py sha256 <64 hex digits>": the hash is replaced. The archived logs record the sha256 of the
   scripts as originally run; the shipped scripts differ from those only in file paths, comments, the library's file name, the
   removal of the hash check mentioned above and the one label;
2. run times: a number followed by "s", after "(" or a space and before ")", ";" or the end of the line (for example "(14s)",
   " 41s;", " 4s"), becomes "<N>s".

## Provenance

An independent check made for the paper's revision (October 2026); re-run for release v1.2.8.
