# classical_locus_checks: machine checks of the classical jumping locus

## What the computation shows

The top-degree part of each sheet's I_3 density is a quartic symbol H = x^4 + b x^2 r^2 + c r^4, and a first-derivative
density F of weight W = s + 1 is classically conserved iff H_xr (F_xx - F_rr) - (H_xx - H_rr) F_xr = 0. These scripts check,
in exact arithmetic (sympy and Python fractions; no floating point), the classical results of the jumping-loci paper:
(A4) the symbols b(t), c(t) of the three sheets, computed over the field Q(t) from the repository's certified I_3 densities
(code/mixed_engine.py, cyl_P4) and from the weight-4 profiles (code/cyl_mo2_profiles.py), equal the table of the paper, lie
on the curves c = b/6, c = b^2/(6 - b), c = 1, and give the sheet-by-sheet consequences of Theorem 1 for every odd spin
s <= 61; (A3) the conservation matrix built from the equation above is bidiagonal with the stated entries at W = 4..22, its
maximal minors are the stated products at W = 6..14, and the rank-drop criterion agrees with 1360 direct exact ranks at
W = 4..30 (largest multiplicity 2); (A2) on the line at infinity the drops, from the gcd of all maximal minors, are exactly
the points of case (3) of Theorem 1 at W = 4..30, the line h_1 = 0 has multiplicity 2, and sheet 2 meets the line at
infinity at t = 11/5, a classical jump exactly at spins 9 and 27 in that range; (A5) the ledger: every isolated classical
jumping point of sheets 2 and 3 at odd spins 5-29 (84 points; sheet 1 has none), each verified by the exact classical
multiplicity 2 of the committed symbol at that point, classified as ladder point, recurrence or rule null, with the quantum
outcome and its status label read from the repository's records (no quantum computation is redone here), and the statement
that for every classical jumping point of odd spin <= 61, beta^2 = 1/q with q even holds iff the point is a ladder point t_k.

## Statements of the paper that rest on it

Paper: papers/jumping-loci/main.tex (quoted with the LaTeX markup simplified).

- Sec. 1 (Setting), paragraph "Symbols", after the table of b(t), c(t) and the curves: "identically in t (checked exactly
  over Q(t) against the certified I_3 densities of [VitchevMain] [lab result])". Source: code/run_a4.py, parts (1), (1b), (2);
  logs/run_a4.log.
- Sec. 3 (The two ladders, the lifting rule and the recurrence), the label of Theorem 1: "[hand derivation; certified at
  weights 6--14; checked by machine at weights 4--30, lab result]". The "checked by machine at weights 4--30" part is this
  folder (run_a2.py, run_a3.py, run_a4.py); "certified at weights 6--14" is the registered certificate
  cyl_classical_hierarchy, not this folder.
- App. A (Proof of Theorem 1), paragraph "Machine checks [lab result]": "The bidiagonal form was checked against the matrix
  built from (C) at W = 4,...,22, the minors at W = 6,...,14, and the proposition against 1360 exact ranks at W = 4,...,30
  (largest multiplicity seen: 2); case (3) from the gcd of all maximal minors on the line at infinity at W = 4,...,30; the
  sheet curves against the certified I_3 densities and the consequences above for odd sigma <= 61." Source: run_a3.py
  parts (1), (3), (4); run_a2.py part (2); run_a4.py.
- Sec. 3, paragraph "Even q, sheet by sheet" (classical part): "On sheet 2, beta^2 = 1/q with q odd occurs classically at
  the spins (2i+1)(3q-2); the lowest are t = 2, 3/2, 4/3 at spins 7, 13, 19". Source: run_a4.py part (2). The quantum
  outcomes quoted in that paragraph come from other records.
- Sec. 3, paragraph "The lifting rule": "For every classical jumping point of odd spin sigma <= 61 other than t = 3,
  beta^2 = 1/q with q even holds iff the point is a ladder point t_k [lab result]". Source: run_a5.py (logs/run_a5.log, the
  line "for every classical point at odd s <= 61: beta^2 = 1/(2k) iff t is a ladder point t_k").
- Sec. 3, same paragraph: "Through spin 19 the classical locus has 35 isolated points on sheets 2 and 3 besides t = 3
  (Appendix B): the nine certified jumps, [...] and 26 points where the rule predicts no jump, of which 9 are certified, 7
  certified scan evidence and 10 pinned lab records. No recorded outcome contradicts the rule." and "its null clause, which
  through spin 29 covers 67 cells, tested at the 26 through spin 19". Source: LEDGER.json, LEDGER.md.
- App. B (The classical points through spin 19 and their quantum outcomes): "Every isolated classical jumping point on
  sheets 2 and 3 at odd spins 5--19 (Theorem 1; sheet 1 has none), each verified by the exact classical multiplicity 2 of
  its symbol", the table (its 35 rows equal the ledger rows with spin <= 19 in t, beta^2, (j, j'), kind, outcome and status),
  the totals, and "At spins 21--29 the same enumeration gives 8 rule-positive cells [...] and 41 nulls". Source: LEDGER.json,
  LEDGER.md.

What the paper takes from elsewhere: every quantum outcome and status label in the ledger is that of the record named in
the row. The ledger was made before two of the spin-21 null cells (t = 11/5 on sheet 3, t = 2 on sheet 2) were computed;
it lists them as untested, and the paper's statement that they were tested rests on a separate computation. For the nine
certified nulls at spins 9-13 the ledger names results/cyl_discriminant_w1012.json, results/cyl_discriminant_w14.json and
code/cyl_oddq_controls.py; the paper cites the registered classification cyl_stratification_w1014 for the exact
classification at weights 10-14 (the same jump sets). The other ledger fields (deck_images, t3,
nonclassical_kernels_w22plus) collect records too; the paper's statements on the deck images and on t = 3 cite their
records directly, not this folder.

Open-tool cross-check. The ledger was made when the classification at weights 6-14 was recorded by computations made with
Mathematica (results/cyl_discriminant_w1012.json, results/cyl_discriminant_w14.json, wolfram/lab/cyl_strata_w*.json), and 156
of its entries name one of them. The same classification has since been computed with open-source tools only: the registered
certificate cyl_stratification_w1014 (results/cyl_stratification_w1014.json) at weights 10, 12, 14, and
results/lab/jumping_loci_records/stratification_w6_w8/ at weights 6 and 8. code/open_tool_crosscheck.py checks every one of
those 156 entries against the open-source records (jump exactly at the listed points, with the listed multiplicity; no jump
elsewhere) and fails on any disagreement; logs/open_tool_crosscheck.log records that all 156 agree. So no outcome in the
ledger depends on Mathematica.

## Files

- README.md: this file.
- run.sh: re-runs the four scripts from the release tree into a new output directory.
- compare.sh: compares the new outputs with the shipped records.
- code/jl_lib.py: shared library (the conservation matrix built from the equation, exact ranks, sheet maps, closed forms).
- code/run_a2.py: check A2 (case (3) at infinity, the line h_1 = 0, the sheets against the line at infinity, t = 11/5).
- code/run_a3.py: check A3 (bidiagonal form, maximal minors, the criterion against 1360 exact ranks).
- code/run_a4.py: check A4 (the symbols over Q(t), the curves, the consequences for odd s <= 61, the repository's classical lists).
- code/run_a5.py: check A5 (the ledger; writes LEDGER.json and LEDGER.md; also prints the deck-image and t = 3 tables).
- code/open_tool_crosscheck.py: checks the 156 ledger entries that name a Mathematica record against the open-source records
  (added for release v1.2.8; see "Open-tool cross-check" above).
- code/vendor/mixed_engine.py, invariant_engine.py, cyl_mo2_profiles.py: byte copies of the repository's code/ files of the
  same names; each script compares their sha256 with the repository files before use.
- code/vendor/SOURCE_SHA256.txt: the sha256 of the three vendored files.
- logs/run_a2.log, run_a3.log, run_a4.log, run_a5.log, open_tool_crosscheck.log: the recorded output of each script.
- LEDGER.json: the 84 ledger rows with every record covering each point (all_records), and the deck-image and t = 3 tables.
- LEDGER.md: the same ledger as tables by spin, with summary counts.

Reading the logs: A2-A10 name the items of the check; R1 and R2 are two outside referee reports on an earlier draft of the
paper, whose closed forms of the recursion and statements on the ledger were tested. Mentions of "Table 1's spin-7 row" refer to
that earlier draft, and REPORT.md (named in LEDGER.md) is the check's report, not part of this folder; the t = 3 data are
in LEDGER.json (key t3) and in logs/run_a5.log (section A6).

Repository files the scripts read (all in the public tree): code/mixed_engine.py, code/invariant_engine.py,
code/cyl_mo2_profiles.py; results/lab/classical/rank_drop_locus_w12.txt, rank_drop_locus_w14_18.txt,
rank_drop_locus_w20_22.txt (run_a4.py); and for run_a5.py: code/certifications.json, results/lab/classical/lifting_out.txt,
the run records results/cyl_spin13_prediction_run.json, cyl_spin13_sol3_sporadic_run.json,
cyl_spin15_stratification_run.json, cyl_spin15_ladder_probe_run.json, cyl_spin17_double_prediction_run.json,
cyl_oddq_controls_run.json, cyl_discriminant_w1012_run.json, cyl_discriminant_w14_run.json, the classification records
results/cyl_discriminant_w1012.json and results/cyl_discriminant_w14.json, the stratification outputs
wolfram/lab/cyl_strata_w6.json, _w8, _w10, _w12, _w14 (lab results), the scan summaries
results/lab/modp_scan/w16_p10007_summary.json, w16_p10009, w18_p10007, w18_p10009 and results/lab/spin19/summary.json, and
seven kernel records in results/lab/fast_kernel3/; and for open_tool_crosscheck.py: results/cyl_stratification_w1014.json and
results/lab/jumping_loci_records/stratification_w6_w8/results/S23_sheet{1,2,3}_W{6,8}.json.

## How to re-run

From the repository root (this folder at results/lab/jumping_loci_records/classical_locus_checks/):

    bash results/lab/jumping_loci_records/classical_locus_checks/run.sh /tmp/classical_locus_rerun
    bash results/lab/jumping_loci_records/classical_locus_checks/compare.sh /tmp/classical_locus_rerun

The output directory must not exist yet. Set PYTHON to choose the interpreter (default python3); run.sh calls it as
"$PYTHON -I -B". Recorded environment: Python 3.12.8, sympy 1.14.0, mpmath 1.3.0; python-flint 0.9.0 was installed, so
sympy used its flint ground types; numpy (2.5.2 installed) is not used. Also needed: bash, timeout and cmp.

Expected cost: about 35 seconds in total on one core (run_a2 about 6 s, run_a3 about 15 s, run_a4 about 1 s, run_a5 about
11 s, open_tool_crosscheck under 1 s) and a peak resident memory below 100 MB (measured for run_a2 to run_a4: 77 MB). Each step runs under
ulimit -v 8000000 (kB) and a timeout (900 s for run_a4, 600 s for open_tool_crosscheck, 1800 s for the others).

## What "reproduced" means

run.sh exits 0 and compare.sh prints IDENTICAL for each of logs/run_a2.log, logs/run_a3.log, logs/run_a4.log,
logs/run_a5.log, logs/open_tool_crosscheck.log, LEDGER.json and LEDGER.md, then "ALL IDENTICAL". The comparison is byte for byte, with no mask: the
scripts print no times, host names or absolute paths (run.sh writes its times to run_info.txt in the output directory,
which is not compared).

## Provenance

An independent check made for the paper's revision (October 2026); re-run for release v1.2.8. The open-tool cross-check was
added for release v1.2.8.
