# Provenance of the independent verification engine

The scripts in this directory were written by an **independent
AI-based review** (Claude, Anthropic) of the first versions of the
two accompanying papers (2026-08-18) — the review credited in the
companion paper's Acknowledgments, which found the c_{1,4} = -25/2
identification.  The engine was implemented **from scratch, from
the papers' stated conventions**, sharing no code with this
repository's stack; its results are concordant with the registered
certificates everywhere the two overlap, and extend the
class-by-class tower/KdV coincidence through charge spin 15
(`06_spin14_15_extension.py`).

The files are imported **verbatim** from the package handed over by
the reviewer (SHA-256 of the original archive:
`d4ca1a418dae6baedeaa50ccce4d3210c8d8808caf4a533076c625d1a239df88`).
Only this provenance note was added.  `README.md` in this directory
is the reviewer's own description of the scripts, their
conventions, and per-script runtimes.

## Status

These scripts are **evidence artifacts, not registered
certificates**: they are not run by `code/run_certifications.py`,
carry no pinned markers, and are shipped so that readers can run a
second, independently written engine against the same claims.  The
registered, fail-closed certificates remain those catalogued in
`code/certifications.json`.

## Running

Pure Python 3 + sympy, one script per check, each printing its own
verdicts (see `README.md` for runtimes).  Scripts 07-08 read the
repository's `results/paperclip_P8.json` and `paperclip_P10.json`;
point `REPO_RESULTS` at the `results/` directory of this
repository.
