# integrable-boundary-certificates

Exact-arithmetic engines and machine certificates for the results of

- *Commuting local integrals of motion for free bosonic CFTs and
  integrable boundary states* (E. Vitchev, arXiv:XXXX.XXXXX), and
- its companion on the enhanced c=1 tower (arXiv:XXXX.XXXXX).

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

## Provenance model

Computational reproducibility is public; historical process
provenance remains archived. Commit identifiers inside certified
artifacts are opaque archival references into the private
development history — never mutated here and not resolvable in this
repository by design. Public verification is content-based:
`RELEASE_MANIFEST.json` pins the SHA-256 of every file,
`EXPORT_RECORD.json` records the export (including a small set of
documented local-path redactions with pre/post hashes), and
`code/public_lint.py` re-verifies hashes, artifact cross-agreement,
and the privacy screen from any clean clone.

## License

MIT (see LICENSE). Copyright (c) 2026 Evgueniy Vitchev.
