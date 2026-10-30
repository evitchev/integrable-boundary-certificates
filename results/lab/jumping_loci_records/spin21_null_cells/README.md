# Spin-21 null cells of the lifting rule (weight 22)

## What the computation shows

The lifting rule of the jumping-loci paper (a conjecture) says that a classical jumping point lifts to a quantum one exactly
when beta^2 = (t-1)/(t+1) = 1/q with q even; its null clause says that the other classical points do not lift. This folder tests
the null clause at two cells of spin 21 (weight 22): **S3**, sheet 3 at t = 11/5 (beta^2 = 3/8, (N, s) = (-21/4, -55/4), the
lattice point (j, j') = (3, 3)) and **S2**, sheet 2 at t = 2 (beta^2 = 1/3, (N, s) = (-7, -16), the point (6, 3)). First,
`code/run_g1.py` computes in exact arithmetic the classical multiplicity at weight 22 of each sheet's weight-4 symbol: 2 at S3
and at S2 (both are classical jumping points), and 1 at the two control fibers (sheet 3 at t = 2, sheet 2 at t = 11/5). Then
the exact quantum kernel at weight 22 (the densities whose charges commute with I_3, modulo total derivatives) is computed at
both cells by the modular fast path `fast_kernel3.py`. It takes ranks modulo the 12 largest primes below 2^61 and reconstructs
the kernel vectors by Chinese remaindering and rational reconstruction. Each reconstructed vector is then checked by exact identities in the engine, which
gives the lower bound, and the nullity modulo the reference prime 2^61 - 1 gives the upper bound. Both cells give
"CERTIFIED genuine dim 1 with 4 prime(s)": the kernel has dimension 1, as at a generic point, so neither classical jump lifts.
In the same series of runs the positive control **G2**, the ladder point t_5 = 13/11 on sheet 3, gave "CERTIFIED genuine dim 2
with 5 prime(s)". Its two kernel vectors are equal to those of the repository's earlier record
`results/lab/fast_kernel3/w22_t1311_sheet3.json` (`code/compare_g2.sh`: "WITNESSES EQUAL"), which shows that the computation
detects a jump at this weight. Status: lab result. The kernel records call themselves lab records of the fast path; they are not
recomputed by a registered certificate.

**Prediction recorded before the kernel runs** (28 October 2026). It is taken from the coupling rule recorded in the
repository on 5 September 2026 ("spin 21 (weight 22): sheet 3 at 13/11; sheet 2 nothing"):
- S3: "CERTIFIED genuine dim 1"; S2: "CERTIFIED genuine dim 1".
- G2: "CERTIFIED genuine dim 2", a condition for running S3 and S2. Also expected: 5 primes, and kernel vectors equal to the
  earlier record's.
- The reading was fixed in advance. Dimension 1 means the null clause holds at that cell. Dimension 2 would have required the
  rank of the top-degree parts of the two kernel vectors: rank 2 would contradict the rule, while rank 1 would mean a purely
  sub-leading extra class, which the rule does not constrain. Any other outcome would be inconclusive.
- Outcome: as predicted at all three cells.

**Statements of the paper that rest on this folder** (jumping-loci paper, revision of October 2026):
- Abstract: "its null clauses have been tested through spin 19 and at two cells of spin 21".
- Section 2 (Record through spin 27): "from spin 21 on, only the points predicted by the lifting rule (...) and two of its null
  cells at spin 21 were computed".
- Section 3 (The two ladders, the lifting rule and the recurrence), *Even q, sheet by sheet*: "The next, t=2 at spin 21, does not
  jump either [lab result]".
- Section 3, *The lifting rule*: "These two, at spin 21, were tested against predictions recorded before the computations:
  sheet 3 at t=11/5 (...) and sheet 2 at t=2 (...), each of classical multiplicity 2 at weight 22. At both the quantum kernel
  has dimension 1, so neither classical jump lifts [lab result]. In the same series of runs the positive control, t_5 on
  sheet 3, gave dimension 2 with the same kernel vectors as the existing record."
- Section 7 (Open problems): "Above spin 19 the rule's null clauses are tested only at the two cells of spin 21".
- Appendix B (The classical points through spin 19 and their quantum outcomes): "41 nulls: t=11/5 on sheet 3 and t=2 on sheet 2
  at spin 21 were tested, with no jump (kernel dimension 1) [lab result]".

## Files

- `README.md`: this file.
- `run.sh`: re-runs every computation into a new directory (see below).
- `compare.sh`: compares a re-run with the shipped records and logs (see below).
- `code/run_g1.py`: the classical multiplicities at weight 22 (S3, S2, the two control fibers and G2's fiber), exact.
- `code/jl_lib.py`: the library of the classical checks: the conservation matrix of a quartic symbol and the lattice of eq. (3).
- `code/compare_g2.sh`: compares G2's record with `results/lab/fast_kernel3/w22_t1311_sheet3.json` field by field, and the
  kernel vectors in canonical form (jq -S), with a control that the comparison can fail.
- `code/vendor/lab/fast_kernel3.py`: the exact kernel program (modular ranks, reconstruction, exact verification).
- `code/vendor/code/mixed_engine.py`: byte copy of the repository's `code/mixed_engine.py` (the engine).
- `code/vendor/code/invariant_engine.py`: byte copy of the repository's `code/invariant_engine.py` (imported by the engine).
- `code/SOURCE_SHA256.txt`: the sha256 of every code file as shipped and of the file the original runs used.
- `out/w22_t115_sheet3.json`: the S3 record (genuine_dim 1, 4 primes, 1 kernel vector); 3944525 bytes.
- `out/w22_t2_sheet2.json`: the S2 record (genuine_dim 1, 4 primes, 1 kernel vector); 3617011 bytes.
- `out/w22_t1311_sheet3.json`: the G2 record (genuine_dim 2, 5 primes, 2 kernel vectors); 9850409 bytes.
- `logs/run_g1.log`: the output of `code/run_g1.py`, the classical gate G1 ("G1 PASSED").
- `logs/w22_t1311_sheet3.log`: the output of the G2 kernel run.
- `logs/w22_t115_sheet3.log`: the output of the S3 kernel run.
- `logs/w22_t2_sheet2.log`: the output of the S2 kernel run.
- `logs/compare_g2.log`: the output of `code/compare_g2.sh` ("WITNESSES EQUAL").

A record has the fields tool, tool_sha256 (the sha256 of the program that wrote it), weight, solution (the sheet), N, s,
genuine_dim, kc, primes_used, seconds (run time), status, and witnesses. Each witness is one kernel vector: charge_v maps a
candidate index to its exact coefficient, and d_coefficients_S holds the total-derivative part. The logs hold the programs' own
output from the original runs. The runner's header and exit lines (machine, paths, times) are left out, and the record's path
is shown as `<folder>/out/...` and `<repository>/results/...`. `code/vendor/lab/fast_kernel3.py` is a byte copy of the program
that wrote the records (tool_sha256 07c3e3f2...). `code/jl_lib.py` differs from the file used in the original run in one comment
line and one docstring line; its code is unchanged (the syntax trees are equal once docstrings are removed), and
`logs/run_g1.log` shows the original's sha256 prefix 9ff3e9d0f544908b. Both hashes are in `code/SOURCE_SHA256.txt`.

## How to re-run

From the repository root, with an output directory that does not exist yet:

    PYTHON=python3 results/lab/jumping_loci_records/spin21_null_cells/run.sh /tmp/spin21_rerun
    results/lab/jumping_loci_records/spin21_null_cells/compare.sh /tmp/spin21_rerun

`run.sh NEWDIR [STEP ...]` runs, in this order, the steps g1, g2, compare_g2, s3 and s2, or only the steps named. It never writes
into this folder. Every step runs under `ulimit -v` and `timeout`: 8000000 KiB and 600 s for g1 and compare_g2, and
32000000 KiB and 5400 s for each kernel run, the limits of the original runs. Under a limit of 16000000 KiB the parallel
eliminations failed to allocate. Before each kernel run, `run.sh` stops if MemAvailable is below 80 GiB (`MIN_MEM_GIB` changes
this). It also records the environment in `NEWDIR/run_meta.txt`.

Environment of the original runs (Linux, x86_64):
- Python 3.12.8 and python-flint 0.9.0, as recorded in the original run logs.
- sympy 1.14.0 and mpmath 1.3.0, installed in the same interpreter (with gmpy2 2.3.1, which sympy uses when present).
- numpy is not used by these runs; the interpreter had numpy 2.5.2.
- jq 1.6, GNU bash, coreutils and sed.
- The repository's `requirements.txt` covers sympy and mpmath; python-flint is an extra need of this folder.

Runtime and memory, as measured in the original runs (4 worker processes, on a shared machine):
- g1: under 1 s.
- compare_g2: under 1 min.
- Each kernel run: 1353-1386 s, about 23 min. The steps were the exact build (about 125 s) and the normal forms (about 135 s),
  then the 12 per-prime eliminations in 4 parallel processes (about 950 s) and the exact verification (about 135 s).
- In all, about 70 min.
- Peak memory of each kernel run: about 60 GiB of RAM in all. Each of the 4 eliminations had about 14 GiB resident and an
  address space of up to 15.3 GiB.

## What "reproduced" means

`run.sh NEWDIR` followed by `compare.sh NEWDIR` prints `IDENTICAL` for each of the eight shipped records and logs, and then
`ALL IDENTICAL` (exit status 0). The comparison applies this mask to both sides; nothing else may differ:
- The records are compared as canonical JSON (`jq -S -c`) without the field `seconds`. Their field `tool_sha256` is accepted
  when it is either 07c3e3f227a719703e8a2e6108083871f14f91d997aa4201bcb03f0d30d4846e (the program that wrote the shipped
  records) or the sha256 of the shipped `code/vendor/lab/fast_kernel3.py`.
- The logs are compared line by line, with each duration (a number directly followed by `s`) masked.
- In the line `record written: ...` the directory is masked, and the 16-hex-digit tool hash is accepted if it is the prefix of
  one of the two hashes above.
- In `logs/run_g1.log`, the 16-hex-digit hash of `code/jl_lib.py` is accepted if it is 9ff3e9d0f544908b (the original) or the
  prefix of the shipped file's sha256.
- In `logs/compare_g2.log` the directories in the first two lines are masked. On line 1 (the new record) its sha256 and size
  are masked, along with its `tool_sha256` (accepted as above) and its `seconds`.

The masked copies of both sides are written to `NEWDIR/compare_masked/`.

## Provenance

An independent check made for the paper's revision (October 2026); re-run for release v1.2.8.
