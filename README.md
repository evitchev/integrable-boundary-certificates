# integrable-boundary-certificates

Exact-arithmetic engines and machine certificates for the results of

- *Commuting integrals of motion for free bosonic CFTs revisited:
  certified densities, spectral rigidity, and the solvability of
  the bilinear commutativity systems* (E. Vitchev, 2026), and
- *The enhanced odd-spin tower at c = 1 is quantum KdV at
  c_{1,4} = -25/2: identification, mechanism, and the Virasoro
  vacuum algebra at the self-dual point* (E. Vitchev, 2026).
- *Opers for the cylindrical integrals of motion* (E. Vitchev, 2026);
  its verification kit is `opers-2026-10/` (see below).
- *Jumping loci of the cylindrical commuting families: arithmetic ladders
  on the curve s^2 = (N-25)(N-1)* (E. Vitchev, 2026); its certificates are
  entries of the registry, and the lab records it cites are in
  `results/lab/jumping_loci_records/` (new in v1.2.8; see below).

The directory `super-kdv-2026-10/` (new in v1.2.6) holds an independent computation of the
cylinder vacuum eigenvalues of the N = 1 super-KdV integrals of motion I_1 ... I_9 in the
Neveu-Schwarz and Ramond sectors (see below).

arXiv identifiers will be added upon announcement.

Pure Python; the dependencies are `sympy` and `mpmath` (see `requirements.txt`). `numpy` is optional:
only two checks of the opers kit (`vir1r_*`) use it. One expensive certificate, `cyl_stratification_w1014`
(new in v1.2.8; run by `--full` only), also needs `numpy`, a C compiler and the FLINT library
(see the header of `code/cyl_stratification_w1014.py`).

## Reproduction

    pip install sympy mpmath
    cd code
    python3 run_certifications.py          # fast checks + hash pins
    python3 public_lint.py                 # provenance + privacy lint
    python3 run_certifications.py --full   # + hours-class checks

`clean_venv_gate.sh` (new in v1.2.7) builds a fresh virtual environment with only the declared
requirements and runs the quick certifications, the opers kit's quick suite with its drills and
the super-KdV checks; it exits 0 only if every line passes (`bash clean_venv_gate.sh`).

The Lean certificates (`lean/`) additionally require a local
Lean 4.24 + Mathlib toolchain; their check fails loudly, not
silently, where the toolchain is absent.

## Lab records of the jumping-loci paper (`results/lab/jumping_loci_records/`)

New in v1.2.8: one subfolder per lab result of the paper that is not a row of its Table 1, each with
a README, the scripts, the records, `run.sh` (re-runs the scripts into a new directory) and
`compare.sh` (compares the output with the records). The folder's README lists them, with their
costs; one of them (`spin21_null_cells/`) needs about 60 GiB of memory per kernel run.

## Verification kit for the opers paper (`opers-2026-10/`)

Self-contained: the archived stage directories behind each claim of the
paper, the certified input tables they compare against (with hashes),
the registered predictions, and the Lean proofs (five files since v1.2.4) with a Lake project.

    cd opers-2026-10
    python3 run_all.py --quick     # 212 quick checks, about 2 h; exits non-zero on any failure
    python3 run_all.py --full      # all 235 checks, about 12 h on one core
    python3 run_all.py --drill     # every check re-run with a tamper that must fail

See `opers-2026-10/README.md` for what each check establishes and what
the kit does not establish. Requires `sympy` and `mpmath`.

## N = 1 super-KdV eigenvalues (`super-kdv-2026-10/`)

Code and results for the vacuum eigenvalues of the local integrals of
motion I_1, I_3, I_5, I_7, I_9 of quantum N = 1 super-KdV on the
cylinder, in the NS and Ramond sectors, exact and symbolic in the
momentum and the background charge.

    cd super-kdv-2026-10
    python3 run_checks.py           # about 50 min; must end with ALL PASS
    python3 run_checks.py --drill   # both checks re-run with a tamper that must FAIL

See `super-kdv-2026-10/README.md` for the conventions, the method and
the results files. Requires `sympy` and `mpmath`.

## Independent verification engine

`independent/` contains a second exact free-boson OPE engine and
check scripts written **from scratch, from the papers' stated
conventions**, by the independent review credited in the companion
paper's Acknowledgments — imported verbatim (provenance and the
original archive hash in `independent/PROVENANCE.md`).  Its
results are concordant with the registered certificates and extend
the tower/KdV class-by-class coincidence through charge spin 15.
These are evidence artifacts, not registered certificates: they
are not run by `run_certifications.py`.

## Verification triangle (Mathematica)

`wolfram/triangle_targets.wls` is a second independent replication,
written by a clean-context implementer from the papers and a hand-off
spec only (no access to this repository's engines, stored results, or
code; the spec states the expected values, so this is an independent
IMPLEMENTATION checked against specified regression targets, not a
blind prediction), on the Wolfram Engine with K. Thielemans'
OPEdefs 3.1 (`wolfram/OPEdefs/`, redistributed unmodified under its
own terms -- see the NOTICE there; not covered by this repository's
MIT license).  It reproduces the headline targets -- the two
Theorem-1 operator identities, the odd-sector commutator constant,
the spin-3 screening locus, the thirteen L(1,0) dimensions, and the
N = 28 cylindrical statements -- and its output
(`results/wolfram_triangle.json`) is compared fail-closed with
`code/wolfram_targets.json` by the EXPENSIVE entry `wolfram_triangle`
(`code/wolfram_triangle.py`), which requires a local `wolframscript`
and fails loudly, never silently, where it is absent.  The
certificate authority stays in the free toolchain.

## Provenance model

Computational reproducibility is public; historical process
provenance remains archived. Commit identifiers inside certified
artifacts are opaque archival references into the private
development history — never mutated here and not resolvable in this
repository by design. Public verification is content-based:
`RELEASE_MANIFEST.json` pins the SHA-256 of every file in the
release (everything except the manifest itself),
`EXPORT_RECORD.json` records the export (documented redactions as
rule labels, counts, and post-redaction hashes), and
`code/public_lint.py` re-verifies hashes, release completeness,
artifact cross-agreement, and the privacy screen from any clean
clone.

## License

MIT (see LICENSE). Copyright (c) 2026 Evgueniy Vitchev.  Exception:
`wolfram/OPEdefs/` is K. Thielemans' OPEdefs, redistributed unmodified
under the terms in its own header (see `wolfram/OPEdefs/NOTICE`).
