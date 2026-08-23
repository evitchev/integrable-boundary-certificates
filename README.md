# integrable-boundary-certificates

Exact-arithmetic engines and machine certificates for the results of

- *Commuting integrals of motion for free bosonic CFTs revisited:
  certified densities, spectral rigidity, and the solvability of
  the bilinear commutativity systems* (E. Vitchev, 2026), and
- *The enhanced odd-spin tower at c = 1 is quantum KdV at
  c_{1,4} = -25/2: identification, mechanism, and the Virasoro
  vacuum algebra at the self-dual point* (E. Vitchev, 2026).

arXiv identifiers will be added upon announcement.

Pure Python; the only dependency is `sympy` (see `requirements.txt`).

## Reproduction

    pip install sympy
    cd code
    python3 run_certifications.py          # fast checks + hash pins
    python3 public_lint.py                 # provenance + privacy lint
    python3 run_certifications.py --full   # + hours-class checks

The Lean certificates (`lean/`) additionally require a local
Lean 4.24 + Mathlib toolchain; their check fails loudly, not
silently, where the toolchain is absent.

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
spec only (no access to this repository's engines, stored results,
or target values), on the Wolfram Engine with K. Thielemans'
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
