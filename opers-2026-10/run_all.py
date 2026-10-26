#!/usr/bin/env python3
"""Verification kit for "Opers for the cylindrical integrals of motion" (E. Vitchev, 2026).

Usage:  python3 run_all.py [--quick | --full] [--drill] [--only NAME[,NAME...]] [--list]

Every path is resolved from this file's location; the working directory does not matter.
Before anything runs, every file of the kit is checked against SHA256SUMS (and the certified input tables
against INPUTS.SHA256).  Each check then runs in a fresh temporary copy of tree/ (so archived logs are never
overwritten), and PASSES only if
  (i) the command exits with the expected code, and
  (ii) its output, with timings removed, equals the archived output of the original run line by line
       (or, where noted, satisfies an explicit numerical criterion).
--drill runs, for every check, a deliberate tamper (a changed input coefficient, a flipped constant, or the
script's own negative-control mode) and PASSES only if every tampered run FAILS the same criterion; where the
tampered run must fail for a stated reason, that reason is required in its output.
Exit status: 0 iff every requested check (and drill) passes; 1 otherwise; 2 on a kit integrity error.
"""
import argparse, hashlib, json, os, re, shutil, subprocess, sys, tempfile, time
from fractions import Fraction

KIT = os.path.dirname(os.path.abspath(__file__))
TREE = os.path.join(KIT, 'tree')
LN = 'results/lab/large_n/'

# ---------------------------------------------------------------- checks
# name: (paper section, kind, cwd relative to tree, argv, output source, reference, expected exit, mode, env)
#   output source 'stdout' or a file the script writes (relative to cwd)
#   reference: archived file relative to cwd, or None with a custom criterion
C = {}
def check(name, section, kind, cwd, argv, ref, out='stdout', code=0, quick=True, env=None, crit=None, post=None, drill=None, note='', mask=None, ref_lines=None, timeout=3600):
    C[name] = dict(timeout=timeout, section=section, kind=kind, cwd=cwd, argv=argv, ref=ref, out=out, code=code, quick=quick,
                   env=env or {}, crit=crit, post=post, drill=drill, note=note, mask=mask or [], ref_lines=ref_lines)

# --- Verification table (Section 6): registered predictions vs certified tables, all five lead comparisons
check('gate', 'Verification (table)', 'registered predictions vs certified tables', LN + 'lead_gates', ['strict_gate.py'],
      'strict_gate_v2.log',
      drill=dict(argv=['strict_gate.py', 'VIR7=results/lab/large_n/lead_gates/tamper_vir7.json'], must=['GATE FAIL', 'spin 11 mismatch 1']),
      note='re-runs the comparators for VIR7 (Sol 1 spins 5-13), VIR4D (Sol 3 spins 11/13), VIR5A (Sol 3 spins 5-13), VIR5C (Sol 3 spin 15, blind), SOL2F (Sol 2 spin 11)')
check('gate_missing_charge', 'Verification (table)', 'gate drill target: a missing prediction must be rejected', LN + 'lead_gates', ['strict_gate.py'],
      'strict_gate_v2.log',
      drill=dict(argv=['strict_gate.py', 'VIR7=results/lab/large_n/lead_gates/tamper_vir7_missing_spin11.json'], must=['GATE FAIL', 'output shape']),
      note='same run as gate; its drill sets the VIR7 t=-2 spin-11 coefficients to null')
# --- Sol 3 operator -> charges, recomputed at t = 2 (k = 2), and the registered VIR4d comparison
check('sol3_operator_t2', 'The opers (Sol 3); Verification', 'recomputation from the operator', LN + 'vir4d', ['d1_predict.py', '2'], 'd1_k2.log',
      post=('same_file', 'prediction_partial_2.json'),
      drill=dict(edit=[('replace', LN + 'vir4d/d1_predict.py', None, None)]),
      note='Gamma-symbol WKB engine; spins 11, 13 at k = 2; output file must be byte-identical to the archived prediction')
check('sol3_vir4d_compare', 'Verification (table)', 'registered predictions vs certified tables', LN + 'vir4d', ['d3_compare.py'], 'd3_compare.log',
      drill=dict(edit=[('bumpat', 'results/lab/anchor11/vev_profile_sol3_w12.json', ['vev', '#1'])], must=[]))
# --- Theorem 1 numerically: three-term ODE == Gamma form, Sol 3 at k = 2 and Sol 2 at k = 10/3, WKB orders 2..8
check('sol3_ode_equiv', 'The opers (Theorem 1)', 'operator identity (formal WKB)', LN + 'vir8', ['a_equiv.py', '3', '2', '8'], 'a_equiv_s3_k2_o8.log',
      drill=dict(argv=['a_equiv.py', '3', '2', '8', 'tamper'], must=['DISAGREE']))
check('sol2_ode_equiv', 'The opers (Theorem 1)', 'operator identity (formal WKB)', LN + 'vir8', ['a_equiv.py', '2', '10/3', '8'], 'a_equiv_s2_k10_3_o8.log',
      drill=dict(argv=['a_equiv.py', '2', '10/3', '8', 'tamper'], must=['DISAGREE']))
# --- Sol 2 registered hold-out
check('sol2_holdout', 'Verification (table)', 'registered prediction vs hold-out table', LN + 'sol2f', ['holdout.py'], 'run_holdout.log', out='run_holdout.log',
      drill=dict(edit=[('bumpat', LN + 'sol2f/inputs_holdout/vev_profile_sol2_w12.json', ['vev', '#1']),
                       ('rehash', LN + 'sol2f/inputs_holdout/vev_profile_sol2_w12.json', LN + 'sol2f/inputs_holdout/ORIGIN_SHA256')], must=[]))
# --- Sol 1: Suzuki + centrifugal operator, zero-parameter check of 882 coefficients
check('sol1_suzuki', 'The opers (Sol 1); Verification', 'operator -> charges vs certified tables (dictionary found post hoc)', LN + 'super3/super3b',
      ['closed_dict.py'], 'closed_dict.log',
      drill=dict(edit=[('bumpat', LN + 'super3/super3b/vev_sol1_w8.json', ['coefficients', '#1'])], must=[]))
# --- Theorem 1 (formal Mellin identities) and Theorem 2 (Ward layer), exact symbolic proofs
check('thm1_mellin', 'The opers (Theorem 1)', 'symbolic identity (proof)', LN + 'proof1', ['p1_equivalence.py'], 'p1_equivalence.log',
      drill=dict(edit=[('replace', LN + 'proof1/p1_equivalence.py', None, None)]))
check('thm2_ward', 'Verification (Theorem 2)', 'symbolic identity (proof)', LN + 'proof1', ['p3_g1_ward.py'], 'p3_g1_ward.log',
      drill=dict(edit=[('bumpat', LN + 'proof1/codex_sol3_projection_loss1.json', ['projection_over_c_a', 'numerator']),
                       ('rehash', LN + 'proof1/codex_sol3_projection_loss1.json', LN + 'proof1/INPUTS.SHA256')], must=[]))
# --- NUM1 (Section 6.1)
check('num1_exact_M1', 'Numerical spectral check', 'numerical (control on an exactly solvable case)', LN + 'num1', ['n1_check.py'], None,
      crit=('num_diffs', 1e-30),
      drill=dict(edit=[('replace', LN + 'num1/n1_check.py', "mp.gamma((2 * l + 3 + e) / 4))", "mp.gamma((2 * l + 3 + e) / 4) * (1 + mp.mpf(10)**-20))")]))
check('num1_doublet_average', 'Numerical spectral check', 'numerical (re-analysis of the archived high-precision determinant data)', LN + 'num1', ['avgpm.py'],
      'run_avgpm.log', out='run_avgpm.log',
      drill=dict(edit=[('bumpat', LN + 'num1/num1_dps60.json', ['data', 'A', 5])], must=[]))
# --- Excited states (Section 8)
check('exc1_level1_I3', 'Excited states (Sol 1)', 'oper vs data (post hoc)', LN + 'exc1', ['o2_compare.py'], 'o2_compare.log',
      drill=dict(argv=['o2_compare.py', 'tamper'], must=['False']))
check('exc1_level1_I5', 'Excited states (Sol 1)', 'registered prediction vs data', LN + 'exc1', ['o3_compare.py'], 'o3_compare.log',
      drill=dict(argv=['o3_compare.py', 'tamper'], must=['FAIL']))
check('exc2_closed_form', 'Excited states (Sol 3)', 'closed form vs data', LN + 'exc2', ['s_closed2.py', 'real'], 's_closed2_real.log',
      drill=dict(argv=['s_closed2.py', 'tamper'], must=['FAIL']))
check('exc2_level1_I5', 'Excited states (Sol 3)', 'registered prediction (blind) vs data', LN + 'exc2', ['u3_compare.py', 'real'], 'u3_compare_real.log',
      drill=dict(argv=['u3_compare.py', 'tamper'], must=['FAIL']))
# --- Boundary (Section 9)
check('br2_commutant', 'Screening systems; Boundary interaction', 'commutant computation', LN + 'br2', ['br2_commutant.py'], 'br2_commutant.log',
      drill=dict(edit=[('replace', LN + 'br2/br2_commutant.py', None, None)]))
check('br2_dims', 'Boundary interaction', 'boundary dimensions', LN + 'br2', ['br2_dims.py'], 'br2_dims.log',
      drill=dict(edit=[('replace', LN + 'br2/br2_dims.py', None, None)]))
# --- Second-order descriptions (Section 5)
check('second1_leading', 'The opers (second-order descriptions)', 'symbolic coefficient', LN + 'second1', ['s2_leading.py'], 's2_leading.log',
      drill=dict(edit=[('replace', LN + 'second1/s2_leading.py', None, None)]))
check('second2_gate1', 'The opers (second-order descriptions)', 'engine control (plant + tamper)', LN + 'second2', ['gate1.py'], 'gate1_plant_tamper.log',
      drill=dict(edit=[('bumpat', LN + 'second2/vev_sol1_w6.json', ['coefficients', '#1'])], must=[]))
check('second2_h2b_t2', 'The opers (second-order descriptions)', 'exact exclusion at the solver output', LN + 'second2_h2b', ['check.py', 'h2b', '2'], 'check_h2b_t2.log', mask=[r'min [-+0-9.e]+ max [-+0-9.e]+'],
      drill=dict(edit=[('bumpat', LN + 'second2_h2b/vev_sol3_w6.json', ['coefficients', '#1'])], must=[]))
check('second2_h2b_plant', 'The opers (second-order descriptions)', 'non-vacuity control (planted dictionary recovered)', LN + 'second2_h2b',
      ['check.py', 'h2b', '2'], 'check_h2b_t2_plant.log', env={'PLANT': '1'}, mask=[r'min [-+0-9.e]+ max [-+0-9.e]+'],
      drill=dict(edit=[('replace', LN + 'second2_h2b/plantcert.py', None, None)], must=[]))
# --- Sol 2 operator at the exact point and further fibres (register item 322); n < 0 engine-zero cells are reported, not passed
check('sol2_ext_exact_point', 'The opers (Sol 2); Verification (table)', 'operator -> charges vs certified tables (exact point t = 1/3)',
      LN + 'sol2_ext_fable', ['ext_fibres.py', '1/3'], 'ext_fibres_run4.log', ref_lines=[0, 1, 28],
      drill=dict(edit=[('replace', LN + 'sol2_ext_fable/ext_fibres.py', None, None)], must=['NO MATCH']))
# --- NUM2 (register item 325): fourth-order Q-functions; spectral Theorem 1 test at t = 2
def num2_fit(tag, quick, section='Numerical spectral check (fourth-order operators)'):
    check(f'num2_fit_{tag.replace("10_3", "t10o3")}', section, 'fit of archived high-precision Q data vs WKB and the registered log E terms',
          LN + 'num2', ['num2_fit.py', tag], f'run_fit_{tag}.log', out=f'run_fit_{tag}.log', quick=quick,
          drill=dict(edit=[('bumpat', LN + f'num2/data_{tag}.json', ['logabsQ', 5, 0], '1e-8')], must=[]))
num2_fit('2_A_60', True); num2_fit('10_3_A_60', True)
for _t in ('2_B_60', '2_C_60', '10_3_B_60', '2_A_120'): num2_fit(_t, False)
check('num2_wrong_ordering', 'Numerical spectral check (fourth-order operators)', 'negative control: the wrong operator ordering must fail (reproduces its archived failing output)',
      LN + 'num2', ['num2_fit.py', '2_A_60_wrong'], 'run_fit_2_A_60_wrong.log', out='run_fit_2_A_60_wrong.log', quick=False,
      drill=dict(edit=[('bumpat', LN + 'num2/data_2_A_60_wrong.json', ['logabsQ', 5, 0], '1e-8')], must=[]))
check('num2_exact_E0', 'Numerical spectral check (fourth-order operators)', 'numerical control: exact E = 0 Q-functions (Meijer G) vs the solver',
      LN + 'num2', ['x1_control.py'], 'run_x1.log', out='run_x1.log',
      drill=dict(edit=[('replace', LN + 'num2/x1_control.py', None, None)], must=['False']))
check('num2_x0_control', 'Numerical spectral check (fourth-order operators)', 'numerical control: generic solver vs the NUM1 determinant', LN + 'num2',
      ['x0_control.py'], 'run_x0.log', out='run_x0.log', quick=False,
      drill=dict(edit=[('bumpat', LN + 'num2/inputs/num1_dps60.json', ['data', 'A', 17], '1e-30'),
                       ('rehash', LN + 'num2/inputs/num1_dps60.json', LN + 'num2/inputs/ORIGIN_SHA256')], must=['False']))
check('num2_t1_zeros', 'Numerical spectral check (spectral Theorem 1 test)', 'numerical: zeros of the projected Q vs the outward-Wronskian determinant', LN + 'num2',
      ['t1_zeros.py'], 'run_t1_zeros.log', out='run_t1_zeros.log', quick=False,
      drill=dict(edit=[('replace', LN + 'num2/t1_zeros.py', None, None)], must=[]))
# --- EXC3 (register item 324): Sol 2 level-one excited states
check('exc3_level1_I5', 'Excited states (Solution 2)', 'registered prediction (blind) vs data', LN + 'exc3', ['u_compare.py', 'real'], 'u_compare_real.log',
      drill=dict(argv=['u_compare.py', 'tamper'], must=['FAIL']))
check('exc3_lead_blind', 'Excited states (Solution 2)', "lead's blind comparison of the registered prediction with its own data blocks", LN + 'lead_exc3_check',
      ['lead_compare_exc3.py'], 'lead_compare_exc3.log',
      drill=dict(argv=['lead_compare_exc3.py', 'tamper'], must=['MISMATCH']))
# --- PROOF2 (v1.2.4): Theorem 2's rank-drop vanishing at losses 2-3 (record side), and the formal-WKB equivalence (Theorem 1) beyond the two archived cases
for _f in ('p1_factor', 'p1_proof', 'p1_proof_b'):
    check(f'proof2_{_f}', 'Spins that are multiples of n', 'symbolic proof (rank-drop vanishing at losses 1-3)', LN + 'proof2', [f'{_f}.py'], f'{_f}.log',
          drill=dict(edit=[('bumpat', LN + 'proof2/codex_sol3_loss2.json', ['projection_over_c_a', 'numerator'])], must=[]))
for _n, _args in (('case1', ['3/2', '1/3', '-', '7']), ('case1_near', ['1500001/1000000', '1/3', '-', '7']), ('case2', ['8/5', '1/3', '-', '7']),
                  ('case3', ['5/3', '1/5', '2/7', '7']), ('case4', ['7/4', '1/3,-2/5', '3/7', '6']), ('case5', ['-6/5', '2/3', '-1/4', '6'])):
    check(f'proof2_general_{_n}', 'The opers (Theorem 1)', 'operator identity (formal WKB), general P and L', LN + 'proof2/vir8copy',
          ['../p2_general.py'] + _args, f'../p2_general_{_n}.log', 
          drill=dict(argv=['../p2_general.py'] + _args + ['tamper'], must=['DISAGREE']))
for _sol, _k, _tag in ((3, '3/2', '3-3_2-8'), (3, '5/3', '3-5_3-8'), (3, '-6/5', '3--6_5-8'), (3, '7/2', '3-7_2-8'),
                       (2, '3/2', '2-3_2-8'), (2, '5/2', '2-5_2-8'), (2, '-6/5', '2--6_5-8')):
    check(f'proof2_aequiv_s{_sol}_k{_k.replace("/", "o").replace("-", "m")}', 'The opers (Theorem 1)', 'operator identity (formal WKB)', LN + 'proof2/vir8copy',
          ['a_equiv.py', str(_sol), _k, '8'], f'../p2_aequiv_{_tag}.log', 
          drill=dict(argv=['a_equiv.py', str(_sol), _k, '8', 'tamper'], must=['DISAGREE']))
# --- PROOF3 (v1.2.4): Theorem 3 (Appendix A) -- the constant omega against the archived Q-function data, and the Hankel-kernel identities
check('thm3_constant', 'Theorem 3 (Appendix A)', 'numerical: (Q4/Q6) m(theta_j)^-1 is one E-independent constant 1/omega (archived NUM2 data, t = 2)',
      LN + 'proof3', ['m3_constant.py'], 'm3_constant.log',
      drill=dict(edit=[('bumpat', LN + 'proof3/t1.json', ['A', '#0', 'Q3', 0], '1e-20')], must=['M3 FAIL']))
check('thm3_kernel', 'Theorem 3 (Appendix A)', 'Hankel kernel, integration by parts, intertwining, frozen zeros, Laplace constant (numerical + symbolic)',
      LN + 'proof3', ['m124_checks.py'], 'm124_checks.log',
      drill=dict(edit=[('replace', LN + 'proof3/m124_checks.py', None, None)], must=[]))
# --- JOINT1 (v1.2.4): the (I_3, I_5) pairing of the level-one excited-state comparisons (mixed trace + explicit pairs)
for _c in range(1, 7):
    check(f'joint1_C{_c}', 'Excited states', 'oper vs data: joint (I_3, I_5) pairing', LN + 'joint1', ['joint1.py', f'C{_c}'], f'C{_c}_real.log',
          drill=dict(argv=['joint1.py', f'C{_c}', 'crossed'], must=['tamper fires']))
# --- AKNS1 (v1.2.4): Solution 1 = quantum AKNS (Fateev-Lukyanov hep-th/0510271)
check('akns_spin3', 'Relation to quantum AKNS', 'FL (95) vs the certified Sol 1 spin-3 table (29 fibres; one fixed dictionary)', LN + 'akns1',
      ['sol1_akns_spin3_check.py', '../../first_order/vev_w4_points.json'], 'sol1_akns_spin3_check.log',
      drill=dict(edit=[('replace', LN + 'akns1/sol1_akns_spin3_check.py', None, None)], must=['0/29']))
check('akns_top', 'Relation to quantum AKNS', 'FL top layer vs the Sol 1 top layer (spins 3-9)', LN + 'akns1', ['sol1_akns_top_check.py'], 'sol1_akns_top_check.log',
      drill=dict(edit=[('replace', LN + 'akns1/sol1_akns_top_check.py', None, None)], must=[]))
# --- NUMTAB (v1.2.4): the spectral evidence table (NUM1-3) with its error estimates
check('numtab', 'Numerical spectral check (evidence table)', 'table of fitted vs predicted coefficients from the archived high-precision data', LN + 'numtab',
      ['numtab.py'], 'run_numtab.log', out='run_numtab.log',
      drill=dict(edit=[('bumpat', LN + 'numtab/inputs/num1/num1_dps60.json', ['data', 'A', 5], '1e-8'),
                       ('rehash', LN + 'numtab/inputs/num1/num1_dps60.json', LN + 'numtab/inputs/ORIGIN_SHA256')], must=[]))
# g1_screenings: check(...) calls for run_all.py (kit v1.2.5), stages VIR1, VIR1-R, VIR2, VIR2-R.
# Paste after the BR2 block of run_all.py; merge DRILLS into drills.json.
# Requires in tree/: results/lab/large_n/{vir1,vir1r,vir2,vir2r}/ (whole stage dirs), code/{mixed_engine,invariant_engine,sparse_linalg}.py,
# results/qt_promotion.json, and the path patches of PATCHES.txt (9 files).
# Paper labels: sec:vir = "The Y-sector is a Virasoro algebra of central charge N" (Proposition 1); sec:screen = "Screening systems";
# sec:boundary = "Boundary interaction".

# --- VIR1 (Fable): Proposition 1, charges in Heis_X (x) Vir_N mod d (formal-N engine vir_lib.py)
check('vir1_engine_selftest', 'The Y-sector is Vir_N (Proposition 1)', 'engine control', LN + 'vir1', ['v0_selftest.py'], 'v0_selftest.log',
      drill=dict(edit=[('replace', LN + 'vir1/vir_lib.py', 'M((), [(2, 2)]): F(-2)}', 'M((), [(2, 2)]): F(-2) + F(1, 1000)}')], must=['FAIL  b', 'SELFTEST FAILED']),
      note='VIR1 C0: formal-N OPE engine vs the record engine code/mixed_engine.py (residues on 189 monomial pairs, d, cyl_P4 at 8 curve points x 3 sols), '
           'Virasoro algebra of T_Y with c = N on 756 triples, loop-independence of the Vir_N states, (T_Y T_Y) composite, dims of Q_4, Q_6; '
           'drill: I_3 transcription with the (22)_Y coefficient -2 + 1/1000 must disagree with cyl_P4 (sec:vir)')
check('vir1_heisvir_I3_I5', 'The Y-sector is Vir_N (Proposition 1)', 'membership test (exact, 120 fibres + special points)', LN + 'vir1', ['v1_main.py'], 'v1_main.log',
      drill=dict(edit=[('replace', LN + 'vir1/vir_lib.py', None, None)], must=['REFUTED']),
      note='VIR1 main (sealed): the W_4 coordinate of I_3 and the three non-Virasoro coordinates of I_5 vanish at 120 generic rational t (all three '
           'solutions), full spin-5 kernels in HV + d at t = -+1/3, 3, 5, -5, 5/3 and N = 1; classical-ladder control (bcd table) reproduced; plants flagged. '
           'Drill (drills.json): Vir_N truncated to words of length <= 2 in vir_lib.py -> I_5 REFUTED (sec:vir)')
check('vir1_heisvir_I5_symbolic', 'The Y-sector is Vir_N (Proposition 1)', 'symbolic in t', LN + 'vir1', ['v4_details.py'], 'v4_details.log',
      drill=dict(edit=[('replace', LN + 'vir1/vir_lib.py', None, None)], must=['DETAILS FAILED']),
      note='VIR1 v4: I_5(t) as a primitive polynomial vector in t (DomainMatrix over Q[t]) for each solution; its three non-Virasoro coordinates are the '
           'zero polynomial ("identically in t" for I_5); limit classes at t = 5/3 and t = 3. Drill: truncated Vir_N gives nonzero coordinates (sec:vir)')
check('vir1_heisvir_I7_kernel', 'The Y-sector is Vir_N (Proposition 1)', 'membership test (own ad I_3 kernel, 16 fibres)', LN + 'vir1', ['v2_spin7.py'], 'v2_spin7.log',
      drill=dict(edit=[('replace', LN + 'vir1/vir_lib.py', None, None)], must=['1 outside']),
      note='VIR1 v2 (post-seal): I_7 from the own ad I_3 kernel mod d at weight 8 lies in HV_8 + d at 16 values of t per solution, incl. kernel-jump points; '
           'plant (44)_Y/1000 flagged. Drill: truncated Vir_N -> I_7 outside (sec:vir)')
check('vir1_heisvir_I7_I11_symbolic_N', 'The Y-sector is Vir_N (Proposition 1)', 'symbolic in N on the certified promoted charges', LN + 'vir1', ['v3_ladder.py'],
      'v3_ladder.log', drill=dict(edit=[('replace', LN + 'vir1/vir_lib.py', None, None)], must=['pointwise outside at 24 of 24']),
      note='VIR1 v3 (post-seal): promoted I_7, I_9, I_11 (results/qt_promotion.json, A(N) + s B(N), an INPUT) have all 10/29/80 non-Virasoro coordinates '
           'zero as polynomials in N, and pointwise at 24 values of t (Sol 2 I_11: 22, two poles of the promoted normalisation); promoted I_7 = own kernel at t = 2. '
           'Drill: truncated Vir_N -> outside at 24 of 24 (sec:vir)')
# --- VIR1-R (cc): explicit-component randomised rank test (independent code; mixed_engine only for the inputs)
check('vir1r_explicit_I3_I5', 'The Y-sector is Vir_N (Proposition 1)', 'explicit components, randomised rank mod p (2 primes x 2 seeds)', LN + 'vir1r', ['vir1r.py'],
      'run_VIR1R.stdout', drill=dict(edit=[('replace', LN + 'vir1r/vir1r.py', None, None)], must=['I_5 (weight 6): OUT']),
      note='VIR1-R: I_3, I_5 IN Heis_X (x) Vir_N mod d in explicit components at N = 25 (t = 0), N = 26 (t = +-1/5), N = 28 control (t = -1/3), all three solutions; '
           'plant I_5 + (33)_Y/1000 OUT. IN is probabilistic (Schwartz-Zippel). Drill: Vir fields truncated to <= 2 factors -> I_5 OUT (sec:vir)')
check('vir1r_explicit_I7', 'The Y-sector is Vir_N (Proposition 1)', 'explicit components, randomised rank mod p', LN + 'vir1r', ['vir1r.py', '--i7'],
      'run_VIR1R_i7.stdout', drill=dict(edit=[('replace', LN + 'vir1r/vir1r.py', None, None)], must=['OUT  votes']),
      note='VIR1-R --i7: I_7 IN at N = 25 (t = 0) only, all three solutions; weight-8 plant (44)_Y/1000 OUT (sec:vir)')
# --- VIR2 (Fable): two-boson screening pairs, Gram data, commutants
check('vir2_gates', 'Screening systems', 'engine gates', LN + 'vir2', ['w0_gates.py'], 'w0_gates.log',
      drill=dict(edit=[('replace', LN + 'vir2/vir2_lib.py', None, None)], must=['FAIL  g1', 'GATES FAILED']),
      note='VIR2 C0: Virasoro algebra of T_phi = (dphi)^2/2 + rho d2phi with c = 1 - 12 rho^2; h(e^{b phi}) = b^2/2 - rho b; (11)(11) - 2(22) -> 4 (T_phi T_phi); '
           'I_3 images even in X. Drill: rho -> 2 rho in T_phi (sec:vir, sec:screen)')
check('vir2_I5_symbolic', 'Screening systems', 'symbolic I_5(t) (input of vir2_pairs_symbolic)', LN + 'vir2', ['w2_i5sym.py'], 'w2_i5sym.log',
      post=('same_file', 'i5_sym.pkl'), drill=dict(edit=[('replace', LN + 'vir2/w2_i5sym.py', None, None)], must=[]),
      note='VIR2 w2: symbolic I_5(t) per solution, verified against the numeric kernel; i5_sym.pkl must be byte-identical to the archived file '
           '(it is the input of w3_symbolic). Drill: Sol 3 I_3 cross coefficient + 1/1000 -> no weight-6 kernel (the script aborts) (sec:screen)')
check('vir2_pairs_symbolic', 'Screening systems; Boundary interaction', 'symbolic in p', LN + 'vir2', ['w3_symbolic.py'], 'w3_symbolic.log',
      drill=dict(edit=[('replace', LN + 'vir2/w3_symbolic.py', None, None)], must=['I_3 conserved identically: False', 'RESULT: FAILED']),
      note='VIR2 w3: S1, S2, S3, S3\' and both Virasoro screenings conserved against I_3 and I_5 identically in p (t = (2+p^2)/(2-p^2)); four negatives '
           'per solution fail. Drill: S1 a^2 = 1 - p^2/5 (sec:screen table; sec:boundary "checked symbolically for I_3, I_5")')
check('vir2_points_gram', 'Screening systems', 'Gram matrices, exponents, known points (symbolic + exact)', LN + 'vir2', ['w5_points.py'], 'w5_points.log',
      drill=dict(edit=[('replace', LN + 'vir2/w4_main.py', None, None)], must=['conserved I_3 False', 'POINTS FAILED']),
      note='VIR2 w5: Gram matrices of {Vir screenings, pairs} symbolic in p (rank 2; sealed Sol 1/Sol 3 Grams reproduced; Sol 3 {v_phi-, S3+-} and '
           '{v_phi+, S3\'+-} from which the D_3^(2)/C_2^(1) Cartan matrices follow by A_ij = 2G_ij/G_ii), linear relations, B:C = a^2:beta^2, and '
           'conservation re-checked at t = 3, infinity, +-5, +-1/3, 5/3, -3. Drill (drills.json, w4_main.py): S1 a^2 = 1 - p^2/5 (sec:screen)')
check('vir2_commutant_w4_w6', 'Screening systems', 'commutant dimension mod d (exact, two fibres)', LN + 'vir2', ['w6_structure.py'], 'w6_structure.log',
      drill=dict(edit=[('replace', LN + 'vir2/w4_main.py', None, None)], must=['commutant dim mod d =  0', 'STRUCTURE FAILED']),
      note='VIR2 w6: commutant mod d of {pair(s), e^{p phi}} among ALL two-boson densities is 1-dim at weights 4 and 6 and contains I_3, I_5, at p = 3/2 '
           '(t = -17) and p = -2/3 (t = 11/7); pair alone 2/5 (Sols 1, 2), 3/7 (Sol 3); e^{p phi} alone 5/13; Sol 3 I_3 in rotated bosons U, W symbolic in p. '
           'Drill: S1 a^2 = 1 - p^2/5 -> dim 0 (sec:screen; sec:boundary "within Heis_X (x) Vir_N the pair singles out the family")')
check('vir2_commutant_w8', 'Screening systems', 'commutant dimension mod d at weight 8 (exact, one fibre)', LN + 'vir2', ['w7_weight8.py'], 'w7_weight8.log',
      drill=dict(edit=[('replace', LN + 'vir2/w4_main.py', None, None)], must=['dim mod d = 0', 'WEIGHT8 FAILED']),
      note='VIR2 w7: weight-8 commutant of {pair(s), e^{p phi}} is 1-dim mod d and equals the promoted I_7, at p = 3/2 (t = -17) only (sec:screen)')
check('vir2_spectrum', 'Screening systems', 'commutant dimensions at weights 2..7 (exact, two fibres)', LN + 'vir2', ['w8_spectrum.py'], 'w8_spectrum.log',
      drill=dict(edit=[('replace', LN + 'vir2/w4_main.py', None, None)], must=['SPECTRUM FAILED']),
      note='VIR2 w8: commutant dims at weights 2..7: Sols 2, 3 [1,0,1,0,1,0]; Sol 1 [1,1,1,1,1,1] with the odd-weight charges X-odd; p = 3/2, -2/3 '
           '(sec:screen "for Solution 1 it also contains X-odd charges at odd weights")')
check('vir2_main', 'Screening systems; Boundary interaction', 'sealed pairs at 30 fibres; I_7 (6 fibres), I_9 (3 fibres) hold-outs; completeness', LN + 'vir2',
      ['w4_main.py'], 'w4_main.log', quick=False, ref_lines=list(range(126)), mask=[r'^total \d+s$'],
      drill=dict(argv=['w4_main.py', 'tamper'], must=['SEALED PREDICTION REFUTED']),
      note='VIR2 main (sealed): S1, S2, S3, S3\' conserved vs I_3, I_5 at 30 fibres (24 real p, 6 imaginary p); hold-outs I_7 at 6 and I_9 at 3 fibres '
           '(promoted charges); joint solution set (Groebner + solve, sympy) = the listed pairs at 6 fibres; controls C1-C5. Drill: the archived scoring '
           'tamper (mirror pairs) (sec:screen "at sample fibres for I_7, I_9"; sec:boundary)')
# --- VIR2-R (cc): independent replication (own formal-N pair engine onwick.py + explicit wick.py)
check('vir2r_engine_control', 'The Y-sector is Vir_N (Proposition 1)', 'engine control (formal N vs explicit components)', LN + 'vir2r', ['ctrl_e.py'], 'run_ctrl_e.log',
      drill=dict(edit=[('replace', LN + 'vir2r/ctrl_e.py', None, None)], must=['T T pole 4: formal-N engine == explicit engine: False']),
      note='VIR2-R control (e): the second formal-N engine equals explicit components at N = 3, 4 (paper: "cross-checked against explicit components at N = 3, 4"). '
           'Drill: loop factor N + 1 in the T T fourth-order pole (sec:vir)')
check('vir2r_replication', 'Screening systems; The Y-sector is Vir_N (Proposition 1)', 'independent replication (exact)', LN + 'vir2r', ['vir2r.py'],
      'run_VIR2R.stdout', drill=dict(edit=[('replace', LN + 'vir2r/vir2r.py', None, None)], must=['I_5 vs pair (beta = -4/5, a^2 = 9/25): V_+ and V_- conserved: False']),
      note='VIR2-R: I_3, I_5, I_7 decompose exactly into Heis_X (x) Vir_N mod d at 9 fibres (3 per solution); each pair conserved vs I_3, I_5, I_7 at the '
           '3 fibres of its solution; commutant of {e^{p phi}, V_+, V_-} 1-dim mod d at weights 4, 6, 8 at the first fibre of each solution (t = -57/7, 33/17, 9/7) '
           'and spanned by I_3, I_5, I_7. Drill (drills.json): Vir words truncated to <= 2 factors -> decomposition fails (sec:vir, sec:screen)')
check('vir2r_posthoc_vir_alone', 'Screening systems; Boundary interaction', 'commutant dimension mod d (post hoc, informational)', LN + 'vir2r', ['posthoc_viralone.py'],
      'run_posthoc_viralone.log', drill=dict(edit=[('replace', LN + 'vir2r/posthoc_viralone.py', None, None)], must=[]),
      note='VIR2-R post hoc, Sol 1 p = 8/5 (t = -57/7): commutant mod d with e^{p phi} alone 5/13/30 and with e^{p phi} + V_+ 2/5/11 at weights 4/6/8 (with the '
           'full pair: 1, in vir2r_replication). The script has no internal assertion; the check is the archived output. Drill: V_+ X-exponent + 1/7 '
           'changes 2 -> 1 (sec:boundary)')
# --- KIT INTEGRITY (v1.2.7, P1b): the logs/SHA256SUMS manifest is part of the kit integrity gate; a stale entry must be caught
check('kit_integrity_manifest', 'Kit integrity (v1.2.7)', 'the integrity gate verifies logs/SHA256SUMS (stale or missing entry must FAIL)', '.', ['kit_integrity_gate.py'],
      None, crit=('integrity',),
      drill=dict(edit=[('stale', 'logs/SHA256SUMS', 'logs/VERSIONS.txt')], must=['KIT INTEGRITY FAIL']),
      note='v1.2.7 (P1b): runs the same verify_integrity() the suite starts with. Drill: one entry of logs/SHA256SUMS set to a wrong hash in a '
           'temporary copy; the gate must then report the mismatch (a missing file or missing manifest fails it too)')
# --- g3_reductions (v1.2.5): the known points (Section "The known points", sec:points) and the Solution 2 n < 0 cells (Verification, exceptional cells)
# Paste after the AKNS1 / NUMTAB block of run_all.py.  DRILLS entries go into drills.json (keys = tree paths, values = [old, new]).
# Two archived scripts carry a documented one-line KIT PATCH (absolute path -> relative to __file__): red1/b_ch.py, nzero1/exact_series_v4.py (PATCHES.txt).
# --- RED1 (anchor t = 3 and N = 0): operator-level limits, exact sympy
check('red1_anchor_operator', 'The known points (anchor t = 3)', 'operator limit (exact symbolic)', LN + 'red1', ['a_op.py'], 'a_op_real.log',
      drill=dict(argv=['a_op.py', 'tamper'], must=['A1 holds: False']),
      note='rescaled Sol 3 operator -> L_inf = oscillator form with the product indicial quadratic; terminating solutions exactly at En = N + 1/2 + rho, '
           'rho in {+-P_X, +-pi}, N = 0..3 (the two c = -2 oscillator ladders), none at generic En; 1F3 vs 1F1 parameter. Drill: drop the 1/2 (the archived tamper mode); A1 must fail. '
           'NOT checked: "not the LCLM" (prose argument in NOTE), the global connection problem (open)')
check('red1_anchor_charges', 'The known points (anchor t = 3)', 'charge-level limit k -> infinity (exact rationals at k = 1e4, 1e8)', LN + 'red1', ['a_ch.py'], 'a_ch.log',
      drill=dict(edit=[('replace', LN + 'red1/a_ch.py', None, None)], must=['NOT as sealed']),
      note='supporting (charge level): oscillator charges = 2^2K B_2K(1/2 + l/2); monic Sol 3 charges -> B(1/2+P_X) + B(1/2+pi), O(1/k). Drill: shift the oscillator centrifugal term by 1/10')
check('red1_n0_operator', 'The known points (N = 0)', 'operator limit (exact operator algebra)', LN + 'red1', ['b_op.py'], 'b_op_real.log',
      drill=dict(edit=[('replace', LN + 'red1/b_op.py', None, None)], must=['NOT as sealed']),
      note='double root of the leading symbol; next coefficient vanishes identically; second-order coefficient (D - d0)^2 - nu^2 with nu^2 = 2(P_X^2 + pi^2) - 1/6 '
           '(= paper eta^2 = 2(P_X^2 + pi^2 - rho^2), rho^2 = 1/12 at N = 0) for Sol 3 (k -> -4) and Sol 2 (k -> -3/2); double-scaling limit of the E term = log Z - e. '
           'Drill: change the symmetric-ordering shift (1+M)/2 = 3 gamma/2 -> gamma in the Sol 3 operator; the next coefficient no longer vanishes. '
           'NOT checked: decoupling of the other branches (not proved; stated in the paper)')
check('red1_n0_charges', 'The known points (N = 0)', 'charge-level limits vs certified tables at t = +-5 (exact rationals, h = 1e-8, 1e-16)', LN + 'red1', ['b_ch.py'], 'b_ch.log',
      drill=dict(edit=[('bumpat', 'results/lab/vev/vev_sol3_w8.json', ['coefficients', '#1'])], must=[]),
      note='supporting (charge level): Sol 3 k -> -4 and Sol 2 k -> -3/2 operator charges -> certified tables at t = 5, -5 and BLZ(M -> 0) at l^2 = 4 Delta, O(h). '
           'Reads results/lab/vev/vev_sol{2,3}_w{6,8,10}.json (already in the kit, INPUTS.SHA256). Drill: bump one Sol 3 spin-7 table coefficient')
# --- RED2 (paperclip t -> infinity, exact point): exact sympy
check('red2_paperclip_sol1_and_exact_local', 'The known points (paperclip; exact point)', 'exact identities', LN + 'red2', ['b1_a1_checks.py'], 'b1_a1_checks.log',
      drill=dict(edit=[('replace', LN + 'red2/b1_a1_checks.py', None, None)], must=['delta): False']),
      note='paperclip: Erdelyi-Kober term identity; the M = 1 determinant variable a = (2l + 3 - E + alpha)/4 is invariant under (alpha, l) -> (alpha - 2 delta, l + delta) '
           '(Sol 1 = oscillator at E - alpha); single member not P_X-even, its even part = doublet average (B_2..B_10). Exact point: local orders of Lambda dw^2 differ '
           '(Suzuki M = -2 vs pillow, pillow orders entered by hand), kappa-odd part covariant (no point map removes it). Drill: a with 2 alpha')
check('red2_paperclip_sol3_tensor', 'The known points (paperclip)', 'tensor product of two Coulomb problems (exact)', LN + 'red2', ['b2_coulomb.py'], 'b2_coulomb.log',
      drill=dict(edit=[('replace', LN + 'red2/b2_coulomb.py', None, None)], must=[]),
      note='builds the tensor product T of the two Coulomb problems (theta^2 - theta - lam_pm) = (x - E x^2); its theta-only part = P4 iff momenta rotated '
           '(no gauge+scale makes T == S: x-parts differ). T is the input that b2_mellin.py copies. Drill: flip the sign of E in one factor')
check('red2_paperclip_sol3_mellin', 'The known points (paperclip)', 'formal Mellin-recurrence identity (exact)', LN + 'red2', ['b2_mellin.py'], 'b2_mellin.log',
      drill=dict(edit=[('replace', LN + 'red2/b2_mellin.py', None, None)], must=[]),
      note='multiplier m(nu) = rho^nu Gamma(nu + 1/2)/Gamma(nu - 1) maps T\'s recurrence to that of Sol 3\'s k = -2 member identically in nu, l0, l1, E (C = 4, mu = rho); '
           'controls: wrong multiplier, unrotated momenta. Contour conditions not proved (as the paper says). Drill: Gamma(nu + 1/3) in place of Gamma(nu + 1/2)')
check('red2_exact_laplace', 'The known points (exact point)', 'exponent obstruction (exact)', LN + 'red2', ['a2_laplace.py'], 'a2_laplace.log',
      drill=dict(edit=[('replace', LN + 'red2/a2_laplace.py', None, None)], must=['lam']),
      note='supporting: Laplace transform L4 of Suzuki(M = -2) has integer exponents {-1,0,1,2} at p = +-i sqrt(E); pillow X-point difference sqrt2 sqrt(-X) is not. '
           'Drill: add -lam chi\'\'/q to L4 (a term of indicial weight); the exponents become lam-dependent, non-integer')
# --- PILLOW1 (exact point): shift-pattern invariants of Mellin recurrences, exact sympy
check('pillow1_mellin_forms', 'The known points (exact point)', 'pillow Mellin forms (exact)', LN + 'pillow1', ['r1_pillow_mellin.py'], 'r1_pillow_mellin.log',
      drill=dict(edit=[('replace', LN + 'pillow1/r1_pillow_mellin.py', None, None)], must=['MINIMAL: (5,']),
      note='the pillow\'s finite recurrences: five terms {0,2,3,4,5} (kappa^2 at 3,5); with the gauge (1 - s^2)^h, h = +-sqrt(-X/2), three terms {0,2,3} (kappa^2 at 3). '
           'Drill: halve the gauge exponent; the three-term form is lost')
check('pillow1_invariants', 'The known points (exact point)', 'shift-ratio invariant (exact)', LN + 'pillow1', ['r1b_invariants.py'], 'r1b_invariants.log',
      drill=dict(edit=[('replace', LN + 'pillow1/r1b_invariants.py', None, None)], must=['invariants found: set()']),
      note='asymmetric gauges give 4-5 terms; the three-term invariant is {-2, 2/3} (pillow) vs {-1/2, 1/3} (Sol 1, t = -1/3) and {-1, 1/2} (Sol 2, t = 1/3). '
           'Imports r1_pillow_mellin (its output is part of the archived log). Drill: gauge exponent sqrt(-X/2) -> (2/3) sqrt(-X/2)')
check('pillow1_doublet', 'The known points (exact point)', 'doublet tensor product shift pattern (exact)', LN + 'pillow1', ['r2_doublet.py'], 'r2_doublet.log',
      drill=dict(edit=[('replace', LN + 'pillow1/r2_doublet.py', None, None)], must=[]),
      note='Suzuki(+alpha) x Suzuki(-alpha) at M = -2: three Mellin terms {8,10,12} (z = 1/x), E at an end; no affine normal form shared with Sym^3(pillow) (8 terms) or '
           'pillow x first order. NOT checked by this script: "depends on alpha only through alpha^2" (paper Sec. 3, l. 178). Drill: E x^1 in place of E x^2 in one member')
check('pillow1_expgauge', 'The known points (exact point)', 'Newton polygon and exponential gauges (exact)', LN + 'pillow1', ['r3_r_expgauge.py'], 'r3_r_expgauge.log',
      drill=dict(edit=[('replace', LN + 'pillow1/r3_r_expgauge.py', None, None)], must=['order >= 1:  True']),
      note='exp(beta s^h), h = 1, 2, 3, on the pillow three-term form: 7, 5, 9 shifts (positive-power exponential gauges); L4 Newton polygon has no slope-0 part. '
           'The Newton points are entered by hand (list "terms"). Drill: give the alpha chi\' term p-weight 0; the regular-type test must flip')
# --- NZERO1 (Solution 2, n < 0 cells t = 4 spin 9, t = 7 spins 3, 9)
check('nzero1_fibres', 'Verification (exceptional cells, Solution 2 at n < 0)', 'operator -> charges at t0 + 1/N vs certified tables (exact rationals)', LN + 'nzero1',
      ['exact_series_v4.py'], 'exact_series_v4.log', quick=False, mask=[r'\[\d+(\.\d+)?s\]', r'Total: \d+s'],
      drill=dict(edit=[('replace', LN + 'nzero1/exact_series_v4.py', None, None)], must=['monic==cert(same t): False']),
      note='at t0 + 1/N, N = 25, 50, 100, 200, the normalised operator coefficient equals the certified charge at that fibre exactly (all three cells), and '
           '|diff(N)|/|diff(2N)| of the largest coefficient of (coefficient - certified(t0)) is ~2 (consistent with 1/N). Removability NOT proved (as the paper says). '
           'Uses the sol2_ext_fable engine (common.py -> cyl_first_order, see REPORT). Drill: momentum scale Cv -> 1.001 Cv (same tamper as the kit\'s sol2_ext drill); an operator tamper k -> 1.001 k also fires but takes ~13 min')
# g2_jumptower: checks for "Spins that are multiples of n" (opers.tex sec:nz, lines ~325-347).
# Paste into run_all.py after the PROOF2 block; merge DRILLS into drills.json.
# Stages: VIR5d (results/lab/large_n/vir5d), VIR8 b-part (large_n/vir8), WUW1, WUW2, and the fusion-law records
# (lab/classical_*fusion*.py + results/lab/fast_kernel{2,3}/ + results/lab/classical/*_out.txt).
_SEC = 'Spins that are multiples of n'
_V5D = LN + 'vir5d'
_TIME_DONE = [r'done in \d+s']

# --- VIR5d Z0(b): loss-1 layer of the operator in closed form; top vanishes <=> loss-1 vanishes <=> 2i <= K-1
check('nz_loss1_identity', _SEC, 'symbolic identity + lattice scan (operator loss-1 layer = Ward loss-1, vanishing on the drop locus)',
      _V5D, ['z0a_f1prime.py'], 'z0a_f1prime.log',
      drill=dict(edit=[('replace', LN + 'vir5d/z0a_f1prime.py', None, None)], must=[]),
      note='VIR5d Z0(b): F1prime == F1 (K = 1..12, five fibres), == the Ward-derived loss-1 amplitudes (K = 2..40, 4095 cells), and '
           'top vanishes <=> loss-1 vanishes <=> 2i <= K-1 on the lattice K = 2..12. Supports "these points coincide with the diagonal rank-drop locus"')

# --- VIR5d Z0(c): the operator's whole un-normalised WKB coefficient vanishes at the drop points (all layers), not at a non-drop lattice point
for _k, _o, _log, _q in (('1/2', '10', 'z0c_k1_2.log', True), ('7', '12', 'z0c_k7.log', True), ('3/2', '12', 'z0c_k3_2.log', True),
                         ('9', '14', 'z0c_k9.log', True), ('5/2', '14', 'z0c_k5_2.log', True), ('1/3', '14', 'z0c_k1_3.log', True),
                         ('5/3', '18', 'z0c_k5_3.log', False), ('-1/2', '8', 'z0c_km1_2.log', True)):
    check(f'nz_op_vanish_k{_k.replace("/", "o").replace("-", "m")}', _SEC, 'operator WKB coefficient, exact (all layers)', _V5D,
          ['z0c_engine.py', _k, _o], _log, quick=_q, mask=_TIME_DONE,
          drill=(dict(edit=[('replace', LN + 'vir5d/z0c_engine.py', None, None)], must=[]) if _k != '-1/2' else
                 dict(argv=['z0c_engine.py', _k, _o, 'stripped'], must=['stripped'])),
          note=f'VIR5d Z0(c): Gamma-symbol operator at k = {_k}, WKB orders <= {_o}; R_2K identically zero exactly at the drop points '
               '(k = -1/2, t = 1/3 is the lattice NON-drop control). Drill: the operator symbol taken at k + 1/1000 (labels unchanged) must not vanish')

check('nz_op_stripped_k1o2', _SEC, 'control: the operator without its Gamma ratio does not vanish at the boundary drop point', _V5D,
      ['z0c_engine.py', '1/2', '10', 'stripped'], 'z0c_k1_2_stripped.log', mask=_TIME_DONE,
      drill=dict(argv=['z0c_engine.py', '1/2', '10'], must=['R identically zero: True']),
      note='VIR5d Z0(d): at k = 1/2 (t = 7/5, i = 2, K = 5 = 2i+1) the stripped symbol (1 - l0^2 T^-2)(1 - l1^2 T^-2) gives R_10 != 0 with empty top: '
           'at the jump-tower points the Gamma (Beta) factor is what makes the coefficient vanish. Drill: the real operator vanishes there')

# --- VIR5d Z1: the certified Solution-3 charge is the limit of the operator's normalised charge at k0 + 1e-8, 1e-16
for _k, _s in (('1', 5), ('3', 7), ('5', 9), ('1/2', 9)):
    for _h in (8, 16):
        _f = f'z1_k{_k.replace("/", "_")}_s{_s}_h{_h}.json'
        check(f'nz_limit_point_k{_k.replace("/", "o")}_s{_s}_h{_h}', _SEC, 'operator charge at k0 + 10^-h (exact rationals)', _V5D,
              ['z1_point.py', _k, str(_s), str(_h)], _f, out=_f,
              drill=dict(edit=[('replace', LN + 'vir5d/z1_point.py', None, None)], must=[]),
              note=f'VIR5d Z1 input: the monic operator charge at k = {_k} + 1e-{_h}, spin {_s}, recomputed; must equal the archived JSON '
                   '(which z1_compare.py reads). Drill: step 2e-h instead of 1e-h must rewrite the file differently')
check('nz_limits_sol3', _SEC, 'directional limits vs certified tables (Sol 3 at t = 5/3, 11/5, 17/7, 7/5)', _V5D, ['z1_compare.py'], 'z1_compare.log',
      drill=dict(edit=[('bumpat', LN + 'vir5d/z1_k5_s9_h8.json', ['coefficients', '#0'])], must=['FAIL']),
      note='VIR5d Z1: max|monic(k0+h) - certified| is O(h) with ratio 1e8 between h = 1e-8 and 1e-16 at all four Sol-3 cells the paper lists '
           '(t = 5/3 spin 5, 11/5 spin 7, 17/7 spin 9, 7/5 spin 9). A limit demonstration in two steps, as the paper says')
# --- VIR8 b: the three-term ODE forms vanish at the same cells; the Sol-2 limit at t = 5/3 (the paper's Solution-2 cell)
check('nz_ode_vanish', _SEC, 'three-term ODE WKB coefficient at the nZ cells (cB = cG = 0)', LN + 'vir8', ['b_nz.py'], 'b_nz.log', code=2,
      drill=dict(edit=[('replace', LN + 'vir8/b_nz.py', None, None)], must=[]),
      note='VIR8 (b): Sol-3 ODE at t = 11/5 (spin 7), 5/3 (5), 7/5 (9) and Sol-2 ODE at 5/3 (5): both master-integral coefficients vanish; exit 2 = '
           'the sealed prediction (cB != 0) was refuted, as archived. Supports "a Beta-function prefactor common to both forms ... has a pole that cancels this zero"')
check('nz_limits_sol2_ode', _SEC, 'directional limits of the ODE families vs certified tables (incl. Sol 2 at t = 5/3)', LN + 'vir8', ['b2_limit.py'], 'b2_limit.log',
      drill=dict(edit=[('bumpat', 'results/lab/vev/vev_sol2_w6.json', ['coefficients', '#0'])], must=[]),
      note='VIR8 b2: Sol-3 ODE -> Sol-3 table at 11/5, 5/3, 7/5 and Sol-2 ODE -> Sol-2 table at 5/3 (spin 5), O(h) at h = 1e-8, 1e-16; '
           'the two cross-checks at 5/3 do not converge (family-dependent limit)')

# --- VIR5d Z2/Z3/Z5: the two-dimensional kernel on the jump tower and the eigenvalue law E_i (i = 1, 2, 3)
for _t, _K, _q in (('5/3', 3, True), ('7/5', 5, True), ('11/5', 4, True), ('17/7', 5, True), ('13/7', 6, True), ('23/9', 6, True), ('9/7', 7, False)):
    _f = f'z2_kernel_t{_t.replace("/", "_")}_K{_K}.json'
    check(f'nz_kernel_t{_t.replace("/", "o")}_K{_K}', _SEC, 'ker ad I_3 mod d at weight 2K (exact) and its vacuum eigenvalues', _V5D,
          ['z2_kernel.py', _t, str(_K)], _f, out=_f, quick=_q,
          drill=dict(edit=[('replace', LN + 'vir5d/z2_kernel.py', None, None)], must=[]),
          note=f'VIR5d Z2: kernel at t = {_t}, spin {2 * _K - 1}; the archived JSON (dimension, classes, eigenvalues) must be reproduced. '
               'Dimension 2 on the jump tower t_i = (2i+3)/(2i+1) (5/3, 7/5, 9/7), 1 at the other drop points')
check('nz_pw_law_i123', _SEC, 'eigenvalue law E_i at i = 1, 2 (read off) and 3 (sealed hold-out), from the kernel JSONs', _V5D, ['z3_pw.py'], 'z3_pw_real.log',
      drill=dict(argv=['z3_pw.py', 'tamper'], must=['P-W (tamper): FAIL']),
      note='VIR5d addendum: antisymmetric class = prod_{j=-i}^{i}(PX+pi-jp)(PX-pi-jp), p^2 = 1/(i+1), all coefficients at i = 1, 2, 3; '
           'no antisymmetric class at 11/5, 17/7, 13/7, 23/9. Drill: p^2 + 1/1000 must fail')
check('nz_kernel_gate_K6_K7', _SEC, 'certified Sol-3 spin-11/13 table charges lie in the computed kernels', _V5D, ['z5_gate_K6_K7.py'], 'z5_gate_K6_K7.log',
      drill=dict(edit=[('bumpat', 'results/lab/anchor13/vev_profile_sol3_w14.json', ['vev', '#1'])], must=[]),
      note='VIR5d gate: the kernels at 13/7, 23/9 (K = 6) and 9/7 (K = 7) contain the certified Sol-3 charge (anchor11/anchor13 profiles)')

# --- Fusion-product law, i = 1..6 (q = 2i+2): fitted at q = 4, 6; predicted and confirmed at q = 8, 10, 12, 14 from the kernel records
_FUS_MASK = [r'VEV computed in \d+s']
check('nz_fusion_fit_q4_q6', _SEC, 'fusion law vs the certified-recipe eigenvalues at q = 4, 6 (i = 1, 2)', 'lab', ['classical_fusion_product.py'],
      '../results/lab/classical/fusion_product_out.txt',
      drill=dict(edit=[('replace', 'lab/classical_fusion_product.py', None, None)], must=['False']),
      note='-prod_{m} Lambda(P + i sqrt2 m beta), Lambda = P^2 - Q^2 + (q-1)^2/(2q), equals the q = 4 (t = 5/3) and q = 6 (t = 7/5) eigenvalues '
           '(hard-coded from results/lab/classical/extra_vev_out.txt), and prints the q = 8 prediction')
_W26_SHA = '67a978bd6565462c8aa0e9bfd9b55dc5e7eca96efbf03145f9e32a483038aa24'   # archived w26_t1513_sheet3.json (60 MB; shipped gzipped)
for _w, _q, _quick in ((14, 8, False), (18, 10, True), (22, 12, True), (26, 14, True)):
    _argv = [f'classical_w{_w}_fusion_test.py'] if _w != 26 else ['kit_gunzip_run.py', '../results/lab/fast_kernel3/w26_t1513_sheet3.json.gz', _W26_SHA, 'classical_w26_fusion_test.py']
    check(f'nz_fusion_q{_q}', _SEC, f'fusion law vs the extra invariant of the weight-{_w} kernel record (i = {_q // 2 - 1})', 'lab',
          _argv, f'../results/lab/classical/w{_w}_fusion_test_out.txt', mask=_FUS_MASK, quick=_quick,
          drill=dict(edit=[('replace', f'lab/classical_w{_w}_fusion_test.py', None, None)], argv=_argv, must=['REFUTED']),
          note=f'q = {_q}, t = {_q + 1}/{_q - 1}, weight {_w}: from the archived witness vectors of the fast-path kernel record (genuine dim 2) the member '
               f'with top (x^2 - r^2)^{_w // 2} is formed, its vacuum eigenvalue computed exactly and compared with the fusion product. Kernel membership '
               'of the witnesses is NOT re-verified (records are lab records pending certification). Drill: Lambda constant +1/(2q) must be REFUTED')

# --- WUW1 (i = 1, t = 5/3) and WUW2 (i = 2, t = 7/5): the second class is W_U W_W at density level
check('nz_wuw1', _SEC, 'density-level identification at i = 1 (sealed tests; archived exit 1: one vacuous control)', LN + 'wuw1', ['wuw1.py'], 'run_wuw1.stdout', code=1,
      drill=dict(edit=[('bumpat', LN + 'wuw1/inputs/vev_sol3_w6.json', ['coefficients', '#0']),
                       ('rehash', LN + 'wuw1/inputs/vev_sol3_w6.json', LN + 'wuw1/inputs/ORIGIN_SHA256')], must=[]),
      note='WUW1: K_6 = ker ad I_3 (weight 6) has dim 2; W_U W_W lies in it and not in span{I_5} + d; eigen(G_6) = -(2/9) E_1. Exit 1 reproduces '
           'the archived failure of the sealed TAMPER 2 control (vacuous for eigenvalues at weight 3, as the paper says: "at i = 1 the eigenvalue match is automatic")')
check('nz_wuw1_uniqueness', _SEC, 'uniqueness of W_U W_W among weight-3 bilinears in K_6 + d (post hoc)', LN + 'wuw1', ['posthoc_products.py'],
      'run_posthoc_products.log', out='run_posthoc_products.log',
      drill=dict(edit=[('replace', LN + 'wuw1/posthoc_products.py', None, None)], must=[]),
      note='WUW1 post hoc: the bilinears sum M_ij f_i(U) g_j(W) in K_6 + d form a 1-dim space, M = (1, 3/2, 1/4)^2 = W_U (x) W_W up to scale. Drill: K_6 without its first basis vector must lose the bilinear')
check('nz_wuw1_vacuity', _SEC, 'why the eigenvalue test is automatic at weight 3 (post hoc)', LN + 'wuw1', ['diag_tamper2.py'],
      'run_diag_tamper2.log', out='run_diag_tamper2.log',
      drill=dict(edit=[('replace', LN + 'wuw1/diag_tamper2.py', None, None)], must=[]),
      note='WUW1 post hoc: eigen(S_3(a d^k U/k!)) = a^3/6 (P^3 - P/4) for every a: the exponent tamper only rescales the eigenvalue')
check('nz_wuw2', _SEC, 'density-level identification at i = 2 (sealed tests)', LN + 'wuw2', ['wuw2.py'], 'run_wuw2.stdout',
      drill=dict(edit=[('replace', LN + 'wuw2/wuw2.py', None, None)], must=[]),
      note='WUW2 (t = 7/5): W^(5) primary of weight 5; among 49 bilinears exactly W^(5) (x) W^(5) lies in K_10 + d; eigen(G_10) = -(81/4800) E_2 '
           '(not automatic at weight 5); primary for T_U + T_W, no representative primary for T_total')
check('nz_wuw2_k10_full', _SEC, 'full K_10 = ker ad I_3 at weight 10 mod d', LN + 'wuw2', ['k10_full.py'], 'run_k10_full.log', out='run_k10_full.log',
      drill=dict(edit=[('bumpat', LN + 'wuw2/inputs/vev_sol3_w10.json', ['coefficients', '#0']),
                       ('rehash', LN + 'wuw2/inputs/vev_sol3_w10.json', LN + 'wuw2/inputs/ORIGIN_SHA256')], must=[]),
      note='WUW2: dim K_10 = 2 (302 = 300 d-exact + 2), certified vev_sol3_w10 at 7/5 in the span, G_10 in K_10 + d')
# --- BR2 (already in the kit tree; not yet a check): the paper's convention check with LVZ's printed paperclip density
check('br2_direct_paperclip', 'Boundary interaction', 'residue identity on a printed density (exact)', LN + 'br2', ['br2_direct.py'], 'br2_direct.log',
      drill=dict(edit=[('replace', LN + 'br2/br2_direct.py', None, None)], must=['False (must be False)', 'DIRECT CHECKS FAIL']),
      note='BR2 (a): LVZ\'s printed paperclip density W4^(sym) (LVZ (71)) is conserved mod d by the four hairpin exponentials (LVZ (57)-(58) and X -> -X) at '
           'n = -18/25, -32/25, -50/169, a tamper of the 6n(n+2) term fails; (b) N = 1 bookkeeping of pair/mirror/Virasoro commutants. '
           'Drill: the (d2Y)^2 coefficient 4n^2+7n+2 -> 4n^2+7n+3 (sec:boundary "as a check of conventions")')
# ===================== g4 (NUM3, EXC4) -- check(...) calls to paste into run_all.py (after the NUMTAB block) =====================
# Merge DRILLS (bottom) into drills.json.  num3_t1_zeros and num3_wkb_pred (full mode) can exceed run_one's default timeout of 3600 s:
# run_all.py should pass a larger per-check timeout for them (e.g. a 'timeout' key in check(), used by run_one) -- see REPORT.
# --- NUM3 (register item 329): third-order Solution 2 operator, Q-function data at t = 2 (k = 2) and t = 9/4 (k = 10/3); Section 6, paragraph
#     "For the third-order equation ..." (sec:num).  Fits are re-run from the archived Q data; the ODE integration is NOT re-run (one point is, as a spot check).
_N3S = 'Numerical spectral check (third-order operator)'
def num3_fit(tag, quick, ext=False):
    script, pre = ('num3_fit_ext.py', 'fitext') if ext else ('num3_fit.py', 'fit')
    nm = f'num3_{pre}_{tag.replace("10_3", "t10o3")}'
    if ext:   # the extended-grid fit uses the data points with E >= 1e4 plus points_ext/: bump an archived points_ext value
        base = tag[:-6] if tag.endswith('_wrong') else tag
        f = LN + f'num3/points_ext/pt_{base}_6' + ('_wrong' if tag.endswith('_wrong') else '') + '.json'
        drill = dict(edit=[('bumpat', f, ['logabsQ', 0], '1e-8')], must=[])
    else:
        drill = dict(edit=[('bumpat', LN + f'num3/data_{tag}.json', ['logabsQ', 5, 0], '1e-8')], must=[])
    check(nm, _N3S, 'fit of archived high-precision Q data vs WKB (sealed bases; ' + ('POST-HOC extended grid E <= 1e9' if ext else 'sealed grid 1e3..1e7') + ')',
          LN + 'num3', [script, tag], f'run_{pre}_{tag}.log', out=f'run_{pre}_{tag}.log', quick=quick, drill=drill,
          note=('the figures the paper quotes (t = 2: 34/25/18/12 digits; t = 9/4: 32/23/16/11, spin 7 to 6; odd K <= 5.1e-8 / 2.9e-5)' if ext else
                'the registered grid: spin 5 1.6-3.5 digits (below 4 everywhere), spin 3 6.0-8.6 (below 7 at t = 9/4 and for one Q_j at t = 2)') + f'; tag {tag}')
for _t in ('2_A_60', '2_B_60', '10_3_A_60', '10_3_B_60'): num3_fit(_t, True, ext=True)      # the paper's figures (extended grid)
for _t in ('2_A_60', '2_B_60', '10_3_A_60', '10_3_B_60', '2_A_120'): num3_fit(_t, True)       # registered grid; 2_A_120 = precision doubling
check('num3_wrong_ordering', _N3S, 'negative control: the ordering x^c T must fail (reproduces its archived failing output, extended grid)',
      LN + 'num3', ['num3_fit_ext.py', '2_A_60_wrong'], 'run_fitext_2_A_60_wrong.log', out='run_fitext_2_A_60_wrong.log', quick=True,
      drill=dict(edit=[('bumpat', LN + 'num3/points_ext/pt_2_A_60_6_wrong.json', ['logabsQ', 0], '1e-8')], must=[]),
      note='paper: "A wrong exponent or the ordering x^c T destroys the agreement" (Sol 2 part); NOTE Controls row "Tamper (b)"')
check('num3_wrong_ordering_sealed', _N3S, 'negative control on the sealed grid (archived failing output)',
      LN + 'num3', ['num3_fit.py', '2_A_60_wrong'], 'run_fit_2_A_60_wrong.log', out='run_fit_2_A_60_wrong.log', quick=True,
      drill=dict(edit=[('bumpat', LN + 'num3/data_2_A_60_wrong.json', ['logabsQ', 5, 0], '1e-8')], must=[]))
check('num3_controls', _N3S, 'numerical controls X0 (NUM1 exact M = 1; archived NUM2 fourth-order point to 5e-59), X1 (exact E = 0 Meijer G, 3rd order, 54.5-56 digits), X1\' (7th order)',
      LN + 'num3', ['controls.py'], 'run_controls.log', out='run_controls.log',
      drill=dict(edit=[('bumpat', LN + 'num3/inputs/data_2_A_60.json', ['logabsQ', 0, 0], '1e-38'),
                       ('rehash', LN + 'num3/inputs/data_2_A_60.json', LN + 'num3/inputs/ORIGIN_SHA256')], must=['False']),
      note='paper: "controls: an archived fourth-order point reproduced to 5e-59, and the exact Meijer-G values at E = 0 to 55-56 digits"')
check('num3_t1_verify', _N3S + ' (spectral Theorem 1 test, Sol 2)', 'verifier of the archived t1_3.json (Q3_j/Q7_j/C_j - 1 recomputed from the archived Q values, C_j exact E = 0) + one T1 energy recomputed from both ODEs',
      LN + 'num3', ['../kit_num3/t1_verify.py'], '../kit_num3/t1_verify.log',
      drill=dict(edit=[('bumpat', LN + 'num3/t1_3.json', ['A', '50.0', 'Q3', 1], '1e-30')], must=['FAIL']),
      note='paper: "Q3_j(E)/Q7_j(-E) = C_j ... for -20 <= E <= 500 at both momentum points and all three j, to at least 36 digits (typically 55)"; recomputes E = 2 at point B (the worst cell, 8.2e-37)')
check('num3_t1_zeros', _N3S + ' (spectral Theorem 1 test, Sol 2)', 'numerical: T1 ratio at 10 energies from both ODEs + zeros of the projected Q vs outward Wronskian (full re-run of t1_zeros3.py)',
      LN + 'num3', ['t1_zeros3.py'], 'run_t1_zeros3.log', out='run_t1_zeros3.log', quick=False,
      drill=dict(edit=[('replace', LN + 'num3/t1_zeros3.py', None, None)], must=[]),
      note='re-runs the sealed T1 script (3534 s under py3.8 on a loaded machine: give it a timeout > 3600 s); the zero at E = -38.52347784... is not quoted in the paper (NOTE Controls)', timeout=10800)
check('num3_spot_point', _N3S, 'spot check: one archived ODE point (t = 2, point A, E = 1e3, dps 60) recomputed from the third-order equation',
      LN + 'num3', ['../kit_num3/spot_point.py', '2', 'A', '60', '0'], '../kit_num3/spot_point_2_A_60_0.log',
      drill=dict(edit=[('bumpat', LN + 'num3/points/pt_2_A_60_0.json', ['logabsQ', 1], '1e-40')], must=['FAIL']),
      note='the data behind every num3_fit* check: this one point is recomputed (12 s) and must equal the archived 60-digit log|Q_j| to 1e-50')
check('num3_assemble', _N3S, 'data provenance: each fitted data_<tag>.json equals the assembly of its 36 archived per-energy point files', LN + 'num3',
      ['../kit_num3/assemble_verify.py'], '../kit_num3/assemble_verify.log',
      drill=dict(edit=[('bumpat', LN + 'num3/points/pt_10_3_B_60_3.json', ['logabsQ', 0], '1e-40')], must=['FAIL']),
      note='ties the data of every num3_fit* check (sealed grid) to the per-energy files, of which num3_spot_point recomputes one')
check('num3_wkb_pred', _N3S, 'WKB charges and the registered log E prediction (R_8 = R_15 = 0) recomputed; wkb_pred3.json must be byte-identical to the registered file',
      LN + 'num3', ['wkb_pred3.py'], 'run_wkb_pred3.log', out='run_wkb_pred3.log', post=('same_file', 'wkb_pred3.json'), quick=False,
      drill=dict(edit=[('replace', LN + 'num3/wkb_pred3.py', None, None)], must=[]),
      note='the WKB targets of every fit and the registered E^-3 log E coefficient zero (paper: "predicted before the data to have coefficient zero"); 2204 s under py3.12, > 3600 s under py3.8 (needs a timeout > 3600 s)', timeout=10800)
for _t in ('2_A_60', '10_3_A_60'):
    check(f'num3_diag_window_{_t.replace("10_3", "t10o3")}', _N3S, 'POST-HOC diagnostic: sealed-basis fit improves by 8-10 digits when E < 1e4 is dropped (the reason for the grid extension)',
          LN + 'num3', ['diag_window.py', _t], 'run_diag_window.log', code=1, ref_lines=[0, 1, 2, 3] if _t == "2_A_60" else [11, 12, 13, 14],
          quick=True,
          drill=dict(edit=[('bumpat', LN + f'num3/data_{_t}.json', ['logabsQ', 20, 0], '1e-8')], must=[]),
          note='archived log holds both tags and a traceback (Emin 1e5 underdetermined; exit 1); the check compares the four fit lines of this tag and the exit code')
check('num3_diag_k8', _N3S, 'POST-HOC: t = 2 fits without the log terms (as good) and the E^-3 collision cell vs lim C_8', LN + 'num3', ['diag_k8.py'], 'run_diag_k8.log',
      out='run_diag_k8.log', quick=True,
      drill=dict(edit=[('bumpat', LN + 'num3/points_ext/pt_2_A_60_6.json', ['logabsQ', 0], '1e-8')], must=[]),
      note='paper: "At t = 2 the non-WKB family meets the WKB series at E^-3 (spin 7), where the coefficient is not a test ... the fits leave it undetermined"')
# --- EXC4 (register item 328): (A) the half-level rule, Delta I_1 = p^2/beta per apparent singularity, confirmed at order 6 (sixth-order form of Sol 3, t = 2);
#     (B) Solution 3 level two: four apparent singularities, five states; I_3 quintic post hoc at t = 2 (1/2,1/3), blind at three registered points.
#     Paper: Section 'Excited states' (sec:exc), paragraph to be added.  The solvers (Singular std/fglm, msolve) are NOT re-run: their archived outputs
#     are verified in pure Python (kit_exc4/*.py), and the polynomial systems they were given are regenerated from the archived scripts.
_EXS = 'Excited states (Solution 3, level two; half-level rule)'
_IBENV = {'IB_LAB': '../../../../lab', 'IB_CODE': '../../../../code'}   # relative to cwd results/lab/large_n/exc4 -> tree/lab, tree/code
for _t, _P, _Q, _blind in (('2', '1/2', '1/3', False), ('9/4', '1/2', '1/3', True), ('2', '1/3', '1/5', True), ('9/4', '2/5', '1/4', True)):
    _tag = 't%s_P%s_Q%s' % (_t.replace('/', '_'), _P.replace('/', '_'), _Q.replace('/', '_'))
    _gtag = 't%s_P%s_Q%s' % (_t.replace('/', '_'), _P.replace('/', '-'), _Q.replace('/', '-'))
    _nm = _gtag.replace('-', '').replace('_', '')
    check(f'exc4_data_level2_{_nm}', _EXS, 'data side: the 6x6 level-two singlet I_3 block (record tool lab/excited_states.py) recomputed; JSON byte-identical',
          LN + 'exc4', ['d2_level2_exc4.py', '3', _t, _P, _Q], f'd2_sol3_{_tag}.log', env=_IBENV, post=('same_file', f'd2_sol3_{_tag}.json'),
          drill=dict(edit=[('replace', LN + 'exc4/vir_lib.py', None, None)], must=[]),
          note=f'charpoly = (level-shift linear factor) x quintic at t = {_t}, (P,Q) = ({_P},{_Q})' + (' (hold-out point; computed after registration)' if _blind else ' (post-hoc point)'))
    check(f'exc4_gb_{_nm}', _EXS, 'pure-Python verification of the archived msolve Groebner basis and Singular fglm eliminants (no solver re-run)',
          LN + 'exc4', ['../kit_exc4/gb_verify.py', _gtag], f'../kit_exc4/gb_verify_{_gtag}.log',
          drill=dict(argv=['../kit_exc4/gb_verify.py', _gtag, 'tamper'], must=['FAIL']),
          note='vdim 40 = 5 configurations x 8 relabellings; I_1 = vacuum + 4; charpoly(M_ee) = quintic^8; quintic == data block'
               + (' == REGISTERED prediction (blind)' if _blind else ' (post hoc)'))
    check(f'exc4_regen_{_nm}', _EXS, 'the four-point polynomial system given to the solvers is regenerated by h_four_q.py (pure Python, before its Singular call); byte-identical',
          LN + 'exc4', ['../kit_exc4/regen.py', f'h_fourq_{_gtag}_four.sing', 'h_four_q.py', _t, _P, _Q, 'four'], f'../kit_exc4/regen_{_gtag}.log', quick=True,
          drill=dict(edit=[('replace', LN + 'exc4/h_four_q.py', None, None)], must=['False']),
          note='ties gb_verify\'s input to the operator class; 45-75 s under both interpreters (one py3.8 run stalled > 12 min and was killed; not reproduced in 10 further py3.8 runs, see REPORT)')
check('exc4_holdout_compare', _EXS, 'registered prediction (blind) vs data: the three hold-out quintics', LN + 'exc4', ['u_compare_exc4.py', 'real'], 'u_compare_exc4_real.log',
      drill=dict(argv=['u_compare_exc4.py', 'tamper'], must=['FAIL']),
      note='PREDICTION_EXC4.json sha256 cff1aead... (asserted against PREDICTION_EXC4.sha256); 0 of 6 coefficients differ at each point; level-shift root present')
check('exc4_lead_blind', _EXS, "lead's blind comparison of the registered prediction with its own level-two blocks (computed first)", LN + 'lead_exc4_check',
      ['lead_compare_exc4.py'], 'lead_compare_exc4.log', drill=dict(argv=['lead_compare_exc4.py', 'tamper'], must=['MISMATCH']),
      note='KIT PATCH (path only): the snapshot PREDICTION_EXC4.json next to the script, same sha256 cff1aeadc77aab1d...')
check('exc4_posthoc_t2', _EXS, 'post hoc: the archived fglm I_3 eliminant at t = 2, (1/2,1/3) equals the quintic factor of the level-two block', LN + 'exc4',
      ['chk_four_gb.py'], 'chk_four_gb.log',
      drill=dict(edit=[('bumpat', LN + 'exc4/d2_sol3_t2_P1_2_Q1_3.json', ['factors', 1])], must=['EQUAL: False']))
check('exc4_pair_regression', _EXS, 'control: the same code with ONE sigma pair gives the EXC2 level-one data quadratic (archived Singular eliminant vs data)', LN + 'exc4',
      ['chk_pair.py'], 'chk_pair.log',
      drill=dict(edit=[('bumpat', LN + 'exc4/d1_sol3_t2.json', ['blocks', 'I3', 'matrix', 0, 0])], must=['equal: False']))
check('exc4_six_onepoint', 'Excited states (half-level rule, order 6)', 'pure-Python verification of the sixth-order one-point solution (archived Singular basis): r41 = -b, Delta I_1 = 1 = half a level',
      LN + 'exc4', ['../kit_exc4/six_verify.py'], '../kit_exc4/six_verify.log', drill=dict(argv=['../kit_exc4/six_verify.py', 'tamper'], must=['FAIL']),
      note='EXC4 A.2 (A1 one point HOLDS): 2 solutions, both r41 = -2/3, I_1 = vacuum + 1, I_3 = 41/64 (not a level-one data root)')
check('exc4_regen_six', 'Excited states (half-level rule, order 6)', 'the sixth-order one-point system is regenerated by f_six_sol3.py (before its Singular call); byte-identical',
      LN + 'exc4', ['../kit_exc4/regen.py', 'f_six_t2_1-2_1-3.sing', 'f_six_sol3.py', '2', '1/2', '1/3'], '../kit_exc4/regen_six.log',
      drill=dict(edit=[('replace', LN + 'exc4/f_six_sol3.py', None, None)], must=['False']),
      note='includes T0: the vacuum Q6 through the Gamma-form engine reproduces the certified Sol 3 I_3 at t = 2')
check('exc4_regen_six_charges', 'Excited states (half-level rule, order 6)', 'the I_1 / I_3 relations of the one-point solutions are regenerated by f_six_charges.py (Gamma-form engine; before its Singular call)',
      LN + 'exc4', ['../kit_exc4/regen.py', 'f_six_charges_1-2_1-3.sing', 'f_six_charges.py', '1/2', '1/3'], '../kit_exc4/regen_six_charges.log',
      drill=dict(edit=[('replace', LN + 'exc4/f_six_charges.py', None, None)], must=['False']),
      note='monic I_1 of the deformed sextic = vacuum - 3 r41/2 (so r41 = -b = -2/3 gives +1)')
# drills.json entries (file -> [old, new], first occurrence replaced) for the ('replace', file, None, None) drills above
# --- full mode only
check('sol3_vir5a_compare', 'Verification (table)', 'registered predictions vs certified tables (all VIR5a cells)', LN + 'vir5', ['e3_compare.py'],
      'e3_compare.log', quick=False, drill=dict(edit=[('bumpat', 'results/lab/vev/vev_sol3_w10.json', ['coefficients', '#1'])], must=[]))
check('second2_h2b_t9o4', 'The opers (second-order descriptions)', 'exact exclusion at the solver output', LN + 'second2_h2b', ['check.py', 'h2b', '9/4'],
      'check_h2b_t9o4.log', quick=False, mask=[r'min [-+0-9.e]+ max [-+0-9.e]+'], drill=dict(edit=[('bumpat', LN + 'second2_h2b/vev_sol3_w6.json', ['coefficients', '#1'])], must=[]))
check('sol2_ext_all_fibres', 'The opers (Sol 2); Verification (table)', 'operator -> charges vs certified tables (14 fibres; engine-zero cells reported)',
      LN + 'sol2_ext_fable', ['ext_fibres.py'] + '1/3 2 3/2 4 -2 1/2 -1/3 13/5 11/2 21/11 7 401/100 399/100 701/100'.split(), 'ext_fibres_run4.log', quick=False,
      drill=dict(edit=[('replace', LN + 'sol2_ext_fable/ext_fibres.py', None, None)], must=['NO MATCH']))
check('num1_shoot', 'Numerical spectral check', 'numerical (zeros of D vs shooting)', LN + 'num1', ['n2_shoot.py'], 'run_n2.stdout', quick=False,
      drill=dict(edit=[('replace', LN + 'num1/n2_shoot.py', None, None)], must=['False']))

# ======== g6_b2 (v1.2.6) ========
# g6_b2: checks for the Solution 1 entry "U(1) coset of N=2" of Section "Screening systems" (opers_npb.tex, sec:screen, l. 109 and
# ll. 122-126).  Stages B2 (results/lab/susy4/boundary/b2) and B2b (results/lab/susy4/boundary/b2b).
# Paste into run_all.py (after the last check block); merge DRILLS into drills.json.
# The paper's statements checked:
#   (a) c_{N=2} = 3 - 12/p^2 and the screening (Gram) data of the a = 4 tower's J-neutral screenings = Solution 1's        -> b2_n2frame_*, b2_level1
#   (b) the two towers are the quantum N = 2 KdV hierarchies with a = 4 (Zeitlin's sl(2|1)^(1)) and a = -2 (KL's sl(2|2)^(1)) -> b2_classical_*
#   (c) I_3, I_5 are not the U(1)-vacuum projections of those towers' charges, t = 2, 7/2, 11/3 (+ 13/5 hold-out)            -> b2_towers_*, b2_sol1_*, b2_level2_*, b2_validate2_*
#   (d) I_3, I_5 (and the X-odd I_2) do not commute with those charges, t = 2, 7/2, 13/5                                    -> b2b_*
_B2S = 'Screening systems (Solution 1: U(1) coset of N = 2)'
_B2 = 'results/lab/susy4/boundary/b2'
_B2B = 'results/lab/susy4/boundary/b2b'
_FRESH = '../../kit_b2/fresh.py'      # kit helper: deletes the archived copy of the script's JSON output, then runs the script unchanged
_DONE = [r'done( \w+)? in \d+(\.\d+)?s']
_JSON_TIMING = [r'"secs": [0-9.]+', r'"rss_MB": \d+']
_FIB = {'2': ('2-3', '2'), '7/2': ('10-9', '7-2'), '11/3': ('8-7', '11-3'), '13/5': ('8-9', '13-5')}   # t -> (s tag, t tag)

# --- (a) N = 2 frame (V1) and Level I
check('b2_n2frame_v1', _B2S, 'N = 2 superconformal OPEs of the two-superfield frame, symbolic in s = p^2 (exact)', _B2, [_FRESH, '--same-json', 'n2frame_v1.json', '--', 'n2frame.py', 'v1'],
      'logs/n2frame_v1.log', mask=_DONE,
      drill=dict(edit=[('replace', _B2 + '/n2frame.py', None, None)], must=['V1 TT/4 = c/2: FAIL', 'FAILURES']),
      note='B2 V1: T, G+, G-, J of the frame eta = (s, -s) with background charge fixed by h(S1) = h(S2) = 1 close into the N = 2 SCA with '
           'c = 3 - 12/s (36 OPE checks, poles 1-5, symbolic in s). Supports "c_{N=2} = 3 - 12/p^2". '
           'n2frame_v1.json rewritten equal as a JSON object. Drill: target c = 3 - 12/s -> 3 - 11/s; the TT and JJ checks must FAIL')
check('b2_n2frame_vertices', _B2S, 'commutation of T, G, g\', J with the five screening vertices at six values of s (exact)', _B2, [_FRESH, '--same-json', 'n2frame_vertices.json', '--', 'n2frame.py', 'vertices'],
      'logs/n2frame_vertices.log', mask=_DONE,
      drill=dict(edit=[('replace', _B2 + '/n2frame.py', 'lam0 = -1 / sv if not isinstance(sv, F) else F(-1) / sv',
                        'lam0 = -2 / sv if not isinstance(sv, F) else F(-2) / sv')], must=['J: S1=n S2=n S3=n S4=n']),
      note='B2: T, G, g\' commute with S1..S4 and W0; J with S1, S2, W0 only (S3, S4 carry J-charge -+2). Supports "U(1) coset": the J-neutral '
           'screenings are S1, S2, W0. n2frame_vertices.json rewritten equal as a JSON object. Drill: background charge -1/s -> -2/s; J must '
           'then fail to commute with S1, S2')
check('b2_level1', _B2S, 'Level I: Gram matrix of the J-neutral screenings = Solution 1\'s, symbolic in s (exact)', _B2, ['level1.py'], 'logs/level1.log',
      drill=dict(edit=[('replace', _B2 + '/level1.py', None, None)], must=['EQUAL: False']),
      note='B2 Level I: after bosonisation and rotation (J = const dH~), Gram(S1, S2, W0) = Gram(V+, V-, e^{p phi}) of Solution 1 '
           '(beta = -p/2, kappa^2 = 1 - p^2/4) identically in s; isometry u = -p phi, v = -2 kappa X. Supports "refers to the central charge '
           'and the screening (Gram) data". Drill: kappa^2 = 1 - s/4 -> 1 - s/3 must give EQUAL: False')
# --- (b) which hierarchies: the classical limits fit the N = 2 KdV flow with a = 4 (Z) and a = -2 (K)
for _tw, _a, _oa in (('Z', '4', '-2'), ('K', '-2', '4')):
    _ot = 'K' if _tw == 'Z' else 'Z'
    check(f'b2_classical_{_tw}', _B2S, f'V2: classical limit of tower {_tw} fits the N = 2 KdV flow with a = {_a} (exact)', _B2, [_FRESH, '--same-json', f'classical_{_tw}.json', '--', 'classical.py', _tw],
          f'logs/classical_{_tw}.log', mask=_DONE,
          drill=dict(argv=[_FRESH, '--same-json', f'classical_{_tw}.json', '--', 'classical.py', _ot], must=[f'a: {_oa}']),
          note=f'B2 V2: classical commutant dims = quantum ones (w2 = 2..8); the quantum weight-4 density at s = 1/j (j = 2..31) reconstructed '
               f'and evaluated at s = 0 equals the classical one; the N = 2 superfield flow fits Mao-Liu (1) with a = {_a} only (unique). '
               f'Supports the labels a = 4 / a = -2 (LM-N2KdV labelling). Drill: the other tower must fit a = {_oa} (the fit discriminates)')

# --- (c) inputs: tower charges (commutants mod d), recomputed at the four fibres; JSON compared with timings masked
for _t, (_st, _tt) in _FIB.items():
    for _tw in ('Z', 'K'):
        for _w, _wa in (('w2_3_4_5_6_7_8', '2,3,4,5,6,7,8'), ('w10_12', '10,12')):
            _f = f'towers_{_tw}_s{_st}_{_w}.json'
            check(f'b2_towers_{_tw}_t{_tt.replace("-", "o")}_{_w.split("_")[0]}{"_8" if _w.startswith("w2") else "_12"}', _B2S,
                  f'tower {_tw} commutant mod d at w2 = {_wa} (t = {_t}, s = {_st.replace("-", "/")}; exact)', _B2,
                  [_FRESH, '--rm', _f, '--', 'towers.py', _tw, _st.replace('-', '/'), _wa], _f, out=_f, mask=_JSON_TIMING,
                  drill=dict(edit=[('replace', _B2 + '/towers.py', '    perturb = None', '    perturb = (2 if tower == "Z" else 3, F(1, 7))')], must=['+(1/7,0)']),
                  note=f'B2 towers: the charges of tower {_tw} (dims and full kernels) that level2.py / b2b.py read; the archived JSON is deleted '
                       'before the run (kit_b2/fresh.py) and must be rewritten identically (secs, rss_MB masked). Drill: the script\'s own '
                       'negative control (last screening momentum + 1/7 along X_0) switched on by a one-line edit; the kernels must change')

# --- (c) Solution 1's certified charges in the coset frame, at the four fibres
for _t, (_st, _tt) in _FIB.items():
    _f = f'sol1_t{_tt}.json'
    _dr = (dict(edit=[('bumpat', 'results/lab/vev/vev_sol1_w6.json', ['coefficients', '#1'])], must=['kit table at t = ' + _t + ' (article P = i P_ours): False'])
           if _t != '11/3' else
           dict(edit=[('replace', _B2 + '/sol1.py', 'E3_repo = sp.expand(CV.vev_density(I3m, sp.Rational(N), P2, Q2))',
                        'E3_repo = sp.expand(CV.vev_density(I3m, sp.Rational(N) + sp.Rational(1, 1000), P2, Q2))')], must=['dictionary: None']))
    check(f'b2_sol1_t{_tt.replace("-", "o")}', _B2S, f'certified Solution 1 I_3, I_5 in the coset frame (t = {_t}; exact)', _B2,
          [_FRESH, '--same-json', _f, '--', 'sol1.py', _t], f'logs/sol1_t{_tt}.log', drill=_dr,
          note=f'B2: certified I_3 (code/mixed_engine.cyl_P4) -> Vir_N -> (X, phi) -> (u, v) by the Level-I isometry: rational and in the '
               f'commutant of Solution 1\'s three (u, v) vectors; momentum dictionary pinned by the exact certified I_3 eigenvalue; I_5 = '
               f'ker ad I_3 at weight 6, its eigenvalue vs the kit table results/lab/vev/vev_sol1_w6.json' +
               (' (kit normalisation singular at t = 11/3: comparison skipped, as archived; Q^6 normalisation used)' if _t == '11/3' else '') +
               f'. sol1_t{_tt}.json (read by level2.py, b2b.py) is deleted and must be rewritten equal to the archived one as a JSON object '
               '(kit_b2/fresh.py --same-json; exit 3 otherwise' + (': at this fibre the archived file has one key in a different position, see REPORT' if _t in ('2', '7/2') else '') + '). Drill: ' +
               ('one coefficient of the kit table bumped by 1/1000 must break the I_5 match' if _t != '11/3' else
                'the certified I_3 eigenvalue evaluated at N + 1/1000 must leave no momentum dictionary (dictionary: None)'))

# --- (c) Level II: projections of the tower charges onto the H~ vacuum vs Solution 1's commutant
for _tw, _t in (('Z', '2'), ('Z', '7/2'), ('Z', '11/3'), ('Z', '13/5'), ('K', '2'), ('K', '7/2'), ('K', '11/3'), ('K', '13/5')):
    _st, _tt = _FIB[_t]
    _f = f'level2_{_tw}_t{_tt}.json'
    check(f'b2_level2_{_tw}_t{_tt.replace("-", "o")}', _B2S + (' (hold-out fibre)' if _t == '13/5' else ''),
          f'Level II: U(1)-vacuum projection of tower {_tw} charges vs Solution 1\'s commutant (t = {_t}; exact)', _B2,
          [_FRESH, '--rm', _f, '--', 'level2.py', _tw, _t], f'logs/level2_{_tw}_t{_tt}.log', mask=_DONE, post=('same_file', _f),
          drill=dict(edit=[('replace', _B2 + '/level2.py', None, None)], must=['8: (True, 2)']),
          note=f'B2 Level II (t = {_t}' + (', the registered hold-out' if _t == '13/5' else '') + f'): the q = 0 top pieces of the projected '
               'spin-3, 4, 5 charges are NOT in Solution 1\'s commutant mod d (spin 2 is: c_2 = (t-1)/(t+3)); pieces of weight <= 6 and constants '
               'are; controls (wrong dictionary s + 1/7, perturbed W0, perturbed Solution 1 vector) fail as archived. Drill: drop Solution 1\'s '
               'third screening e^{-u} from the comparison; the top pieces then lie in the pair commutant (the failure is due to the third screening)')
# --- v1.2.7 (P2): the B2 row labelled "(must be True)" is an informational probe, not a control
for _tw, _t in (('Z', '2'), ('K', '2')):
    _st, _tt = _FIB[_t]   # the expected output is the shipped logs/level2_w8_agreement_{Z,K}.log (covered by SHA256SUMS)
    check(f'b2_level2_w8_agreement_{_tw.lower()}', _B2S + ' (v1.2.7 relabel)',
          f"the w8 real-projection row of level2 (tower {_tw}, t = {_t}) reads False and agrees with the same log's w2 = 8 membership", _B2,
          ['../kit_b2/level2_w8_agreement.py', _tw], f'logs/level2_w8_agreement_{_tw}.log', mask=_DONE,
          drill=dict(edit=[('delete', f'{_B2}/logs/level2_{_tw}_t{_tt}.log')], must=['AGREEMENT CHECK FAILED']),
          note=f'B2 Level II (v1.2.7, P2): the row "real projection at s in the (X, phi) frame: D_0 in span(I_3^cert) + Im d" printed by level2.py is the '
               'real-dictionary projection at w2 = 8, printed for consistency with the rational membership test -- it is not a positive control. This check '
               "enforces exactly that agreement on the archived run: the row reads False AND the same log's w2 = 8 membership line ends ALL: False. "
               'Drill: the archived level2 log removed -> the check cannot establish the agreement and must fail.')
for _t, (_st, _tt) in _FIB.items():
    check(f'b2_validate2_t{_tt.replace("-", "o")}', _B2S, f'end-to-end validation of the projection machinery (t = {_t}; exact)', _B2,
          ['validate2.py', _t], f'logs/validate2_t{_tt}.log',
          drill=dict(edit=[('replace', _B2 + '/coset.py', None, None)], must=['False   diff =']),
          note='B2 validate2: (1) bosonisation of neutral fermion words preserves NS vacuum eigenvalues; (2) for every charge of Z (w2 = 4..12) and '
               'K (4, 8, 12) the eigenvalue of the original two-superfield density equals that of its projection at q = -2P_1/s. '
               'Drill: H~ metric (s-4)/s -> (s-3)/s in coset.py must break (2)')

# --- (d) B2b: commutation of Solution 1's (lifted) charges with the towers' charges
for _tw in ('Z', 'K'):
    for _t in ('2', '7/2', '13/5'):
        _st, _tt = _FIB[_t]
        _f = f'b2b_{_tw}_t{_tt}.json'
        check(f'b2b_{_tw}_t{_tt.replace("-", "o")}', _B2S, f'B2b: [oint L_s, oint I_r] mod d, tower {_tw}, t = {_t} (exact)', _B2B,
              [_FRESH, '--same-json', _f, '--ignore-key', 'secs', '--', 'b2b.py', _tw, _t], f'logs/b2b_{_tw}_t{_tt}.log', mask=_DONE,
              drill=dict(edit=[('bumpat', _B2B + f'/inputs/sol1_t{_tt}.json', ['I3_uv', 0, 1])], must=['[L_2, L_3] = 0 mod d: False']),
              note=f'B2b: residues of L_s (Solution 1\'s I_2, I_3, I_5 lifted H~-free to (u, v, H~)) against tower {_tw}\'s I_r (r = 2, 3, 5' +
                   (', I_2 of Z only' if _tw == 'Z' else '') + ') are NOT total derivatives for every pair; consistency rows (L\'s among '
                   'themselves; I\'s among themselves in the fermionic and the bosonised frame) commute; kappa + 1/7 control. The archived b2b JSON is '
                   'deleted and must be rewritten equal as a JSON object (timings "secs" ignored). Inputs: b2b/inputs/ '
                   '(byte copies of b2/ files, hashes in b2b/inputs/INPUTS.SHA256 and b2b/SHA256SUMS). Drill: one coefficient of the certified '
                   'I_3 in (u, v) bumped by 1/1000 must break [L_2, L_3] = 0')

# ======== g5_spec1 (v1.2.6) ========
# g5_spec1: check(...) calls for run_all.py (kit v1.2.6), stage SPEC1 (register items 387-388).
# Paste after the PROOF3 block (thm3_constant, thm3_kernel) of run_all.py.  No drills.json entries: every drill edit carries its old/new
# strings inline (('replace', file, old, new)) or is a bumpat/rehash.
# Requires in tree/ (FILES.txt): results/lab/large_n/spec1/{code/spec1_{lib,s0,s0b,s1,s2,s2c_gauss,s3,s4_checks,ddt}.py,
# deps_snapshot/{ode,ops}.py, DEPS.SHA256, inputs/proof3/t1.json, out/s3/*.json} and the kit helper dir results/lab/large_n/kit_spec1/
# (verify_s3.py, spot_s3.py, ORIGIN_SHA256, their reference logs, and ref/: the archived SPEC1 logs with the ARGV line removed and the
# seat's absolute workspace prefix replaced by "<stage>/", see README_SECTION.md).  No archived file is patched.
# Paper: papers/opers-2026-10/npb/opers_npb.tex, Theorem thm:spectral, the paragraphs after it, and Appendix A (app:proof3).
SP1 = LN + 'spec1'
SP1_MASK = [r'^(script|imported) (/|<)\S+', r'^python \S+ sympy \S+ mpmath \S+$', r'time ODE \S+ Gamma \S+', r'^total wall time \d+s$']
SP1_SEC = 'Spectral equivalence (Theorem thm:spectral, App. A)'
_lib = SP1 + '/code/spec1_lib.py'
_g_old = '    return (sp.Rational(pt["s"]) - r) / sp.Rational(pt["sigma"]) + sp.Rational(1 - int(pt["a"]), 2)\n'
_g_new = '    return (sp.Rational(pt["s"]) - r) / sp.Rational(pt["sigma"]) + sp.Rational(1 - int(pt["a"]), 2) + sp.Rational(1, 7)\n'
_w_old = '    return sigma ** (g_of(pt, r) - sp.Rational(1, 2)) * sp.sqrt('
_w_new = '    return sp.Rational(8, 7) * sigma ** (g_of(pt, r) - sp.Rational(1, 2)) * sp.sqrt('
# --- S0: scope and the no-resonance hypothesis (H) at the momentum points of the numerics
check('spec1_scope', SP1_SEC + '; scope', 'exact (rational arithmetic)', SP1, ['code/spec1_s0.py', 'S0'], '../kit_spec1/ref/s0.log', mask=SP1_MASK,
      drill=dict(edit=[('replace', SP1 + '/code/spec1_s0.py', '("B", (F(9, 10), F(1, 5)))', '("B", (F(9, 10), F(9, 10)))')], must=['(H) FAILS']),
      note='SPEC1 S0: Sol 3 at t = (3a+2)/(a+2) (d = 4, r = 0) and Sol 1 at t = 2a+1 (d = 2, r = l0) are in the scope of thm:spectral for a = 1..6 '
           '(deg L = 1, sigma > 0), Sol 2 is not (deg L = 2); (H) holds at both momenta (1/3, 3/5), (9/10, 1/5) at the points t = 2, 5/3, 11/5, 3, 5. '
           'Drill: momentum B with l1 = l0 (coincident exponents) -> (H) FAILS')
check('spec1_scope_sol1_a3', SP1_SEC + '; scope', 'exact (rational arithmetic)', SP1, ['code/spec1_s0b.py', 'S0b'], '../kit_spec1/ref/s0b.log', mask=SP1_MASK,
      drill=dict(edit=[('replace', SP1 + '/code/spec1_s0b.py', '("B", (F(9, 10), F(1, 5)))', '("B", (F(9, 10), F(9, 10)))')], must=['(H) FAILS']),
      note='SPEC1 S0b: (H) at Sol 1 a = 3 (t = 7, N = 1/2) at both momenta (the point missing from the S0 table). Drill: l1 = l0 at B -> (H) FAILS')
# --- S1: Lemma 2 / Proposition A.1 identities (symbolic in a, mu, r, s), frozen lattice, intertwining on monomials, Hankel quadrature
check('spec1_kernel', SP1_SEC, 'symbolic identities (general a) + quadrature', SP1, ['code/spec1_s1.py', 'S1'], '../kit_spec1/ref/s1.log', mask=SP1_MASK,
      drill=dict(edit=[('replace', _lib, _g_old, _g_new)], must=['[FAIL] (d)']),
      note='SPEC1 S1: m(mu+sigma)/m(mu) = s - mu - r - c/2 and m(mu)/m(mu-n) = G_a(s-mu-r) symbolically in a, sigma, s, r, mu (and expanded, a = 1..8); '
           'frozen exponents sigma g + sigma j (j < a) are zeros of m and the frozen Frobenius lattice misses the other frozen exponents (a = 1..8); '
           'L_O K = x^n K x^-n L_Gamma on monomials at six points symbolically in mu, l0, l1, E; K x^mu = m(mu) x^mu and Lemma 2(i) by 50-digit '
           'quadrature at three points; each with its own tamper. Drill: g + 1/7 in spec1_lib.g_of -> the intertwining (d) fails')
# --- S2: beta_O, beta_Gamma, the power and saddle identities, omega(t = 2), and the closed-form omega at E = 0 (Meijer G, 12 points x 2 momenta)
check('spec1_constant', SP1_SEC, 'symbolic (beta, saddle) + numerical E = 0 test of the closed-form omega (dps 60, 120)', SP1, ['code/spec1_s2.py', 'S2'],
      '../kit_spec1/ref/s2.log', mask=SP1_MASK,
      drill=dict(edit=[('replace', _lib, _w_old, _w_new)], must=['[FAIL] S2b', '[FAIL] S2c sealed omega']),
      note='SPEC1 S2: the exponents beta_O, beta_Gamma of thm:spectral by exact conjugation (a = 1..6, both solutions), the power identity '
           '-g q_O + beta_G(1 + q_O/sigma) + q_O/2 = beta_O and the saddle identities (symbolic); omega(t = 2) = 3^(-1/6) pi^(-1/2) exactly; the generic '
           'Meijer-G E = 0 formulas reproduce the archived t = 2 closed forms (deps_snapshot/ops.py); omega Q^O(0) = m Q^Gamma(0) at a = 1..6, both '
           'solutions and momenta. Drill: omega * 8/7 (the paper\'s tamper) in spec1_lib.omega_sealed')
check('spec1_e0_gauss', SP1_SEC + '; E = 0 identity', 'symbolic proof (Gauss multiplication) + numerical ties', SP1, ['code/spec1_s2c_gauss.py', 'GAUSS'],
      '../kit_spec1/ref/s2c_gauss.log', mask=SP1_MASK,
      drill=dict(edit=[('replace', SP1 + '/code/spec1_s2c_gauss.py', '    logw = (g - sp.Rational(1, 2)) * Ls + LD / 2 - L2p / 2 - Ld1 / 2',
                        '    logw = (g - sp.Rational(1, 2)) * Ls + LD / 2 - L2p / 2 + Ld1 / 2')], must=['[FAIL] (G6) log(Q^O_j/Q^G_j)']),
      note='SPEC1 S2c: the exact E = 0 identity of thm:spectral: the a frozen Gamma factors and one ODE factor form a complete Gauss product (G1-G3), '
           'and log(Q^O_j(0)/Q^G_j(0)) - log m(theta_j) + log omega = 0 identically in s, r, theta, sigma, d, a (G6); Gauss multiplication spot-checked '
           '(dps 100), the closed form tied to spec1_lib at four points. Drill: (d-1) moved to the numerator of omega\'s square root -> (G6) residual')
check('spec1_saddle', SP1_SEC + ' (App. A, Prop. A.2)', 'symbolic (exponents, saddle) + numerical (descent constant, Laplace rate)', SP1,
      ['code/spec1_s4_checks.py', 'S4'], '../kit_spec1/ref/s4chk.log', mask=SP1_MASK,
      drill=dict(edit=[('replace', SP1 + '/code/spec1_s4_checks.py', '    qG, qO = c / D, c / (d - 1)\n', '    qG, qO = c / D, c / d\n')],
                 must=['[FAIL] (C1) ODE E-correction power']),
      note='SPEC1 S4 checks: (C1) p = (a+1)/(d+a), the E-correction powers -sigma and sigma(p-1) < 0, q_O = p sigma/(1-p); (C2) tau* and the phase '
           'tau* (v - v^p/p), h(1), h\'(1) = 0, h\'\'(1) = 1-p (symbolic); (C3) the descent function R(y) and a positive lower bound of R/min(y^2, y^p) on a '
           'grid (six p; a grid check, not a proof); (C4) the Laplace main term (2 pi (1-p))^(-1/2) with O(1/lambda) error by quadrature. '
           'Drill: q_O = c/d -> (C1) fails')
check('spec1_ddt', 'Solution 1 (relation to DDT at t = 3); ' + SP1_SEC, 'symbolic (operator identities)', SP1, ['code/spec1_ddt.py', 'DDT'],
      '../kit_spec1/ref/ddt.log', mask=SP1_MASK,
      drill=dict(edit=[('replace', SP1 + '/code/spec1_ddt.py', '    kappa = 4 / (3 * sp.sqrt(3))', '    kappa = 4 / (3 * sp.sqrt(2))')],
                 must=['[FAIL] (i) y^2']),
      note='SPEC1 DDT (information in the stage; supports the paper\'s a = 1 remark): with x = mu y^(2/3), mu^6 = 3, the Sol 1 a = 1 ODE is the DDT sextic '
           'with alpha = -(4/3) l0, l + 1/2 = 2 l1/3, E_S = -kappa E, kappa = 4/(3 sqrt 3); the gauge x^(-l0/3) maps the Gamma form to DDT\'s third-order '
           'equation with exponents g = e + l0/3; composite E_D = kappa^-1 E_S. Checked on x^nu symbolically. Drill: kappa = 4/(3 sqrt 2)')
# --- S3 numerics: positive control, archived-output verifier (quick), spot recomputation (quick), and the five full runs (full mode)
check('spec1_s3_control', SP1_SEC + '; numerics', 'numerical control (solver vs the archived t = 2 Q data)', SP1, ['code/spec1_s3.py', 'S3', 'Sol3', '2', 'ctrl'],
      '../kit_spec1/ref/s3_ctrl.log', mask=SP1_MASK,
      drill=dict(edit=[('bumpat', SP1 + '/inputs/proof3/t1.json', ['A', '#0', 'Q3', 0], '1e-20')], must=['[FAIL] S3 positive control']),
      note='SPEC1 S3 control: the solver path of spec1_s3.py at Sol 3 a = 2 (t = 2) reproduces the archived PROOF3/NUM2 Q3, Q6 (t1.json, '
           'byte-identical to proof3/t1.json) at 10 (E, momentum) pairs to >= 45 digits. Drill: one archived Q3 value * (1 + 1e-20)')
check('spec1_s3_archive', SP1_SEC + '; numerics (five points)', 'verifier of the archived S3 outputs (no ODE solved)', LN + 'kit_spec1', ['verify_s3.py'],
      'verify_s3.log',
      drill=dict(edit=[('bumpat', SP1 + '/out/s3/Sol3_a3.json', ['A', '37', 'QO', 0], '1e-8'),
                       ('rehash', SP1 + '/out/s3/Sol3_a3.json', LN + 'kit_spec1/ORIGIN_SHA256')], must=['VERIFY S3 ARCHIVE: FAILURES']),
      note='KIT HELPER verify_s3.py on the archived out/s3/*.json (sha256 = the stage\'s sealed SHA256SUMS): recomputes omega and m(theta_j) exactly and '
           'the deviations at all five points, 8 energies, 2 momenta (equal to the archived ones), and the paper\'s figures: < 1e-30 (criterion), about '
           '1e-56 or better (worst 1.07e-56), better than 1e-64 except at E = 0 for Sol 1 (there 1.6e-62 .. 1.1e-56), the E = 0 Meijer-G comparison of '
           'the same size, and the omega*8/7 and multiplier-interchange tampers fail at every energy. Not the dps-140 repeat (full checks). '
           'Drill: one archived Q^O value * (1 + 1e-8), hash refreshed')
check('spec1_s3_spot', SP1_SEC + '; numerics (five points)', 'recomputation from the operators at two (E, momentum) pairs per point', LN + 'kit_spec1',
      ['spot_s3.py'], 'spot_s3.log',
      drill=dict(edit=[('replace', SP1 + '/code/spec1_s3.py', '    shift = sp.Integer(0) if wrong_order else c / 2\n', '    shift = sp.Integer(0)\n')],
                 must=['SPOT S3: FAILURES']),
      note='KIT HELPER spot_s3.py: re-solves both operators (archived spec1_s3.qpair, dps 70) at E = 7 (momentum A) and E = 500 (B) at the five points; '
           'the Q\'s equal the archived values to >= 60 digits, the relation holds to < 1e-60, the two tampers fail. Drill: the operator ordering '
           'x^c T instead of x^(c/2) T x^(c/2) in L_O (the stage\'s T-order) -> relation and archive agreement fail')
# timeouts: about 2x the slowest measured run (py3.12.8, pure-Python mpmath: 3272, 3574, 3905, 5653 and ~10000 s; py3.8.10 with gmpy is ~25% faster)
for _sol, _a, _to, _lab in (('Sol3', 3, 21600, 't = 11/5, generic'), ('Sol1', 3, 12000, 't = 7, generic'), ('Sol3', 1, 8000, 't = 5/3'),
                            ('Sol1', 1, 7000, 't = 3'), ('Sol1', 2, 7500, 't = 5')):
    check(f'spec1_s3_{_sol.lower()}_a{_a}', SP1_SEC + '; numerics (five points)', 'numerical (full re-run, dps 70 + dps-140 repeat)', SP1,
          ['code/spec1_s3.py', 'S3', _sol, str(_a)], f'../kit_spec1/ref/s3_{_sol.lower()}_a{_a}.log', mask=SP1_MASK, quick=False, timeout=_to,
          drill=dict(edit=[('replace', _lib, _w_old, _w_new)], must=['[FAIL] S3']),
          note=f'SPEC1 S3 at {_sol} a = {_a} ({_lab}): Q^O_j(E) = (m(theta_j)/omega) Q^Gamma_j(-E) with the closed-form omega at E in {{-20, -5, 0, 1, 7, '
               f'37, 100, 500}}, both momenta, dps 70 (< 1e-30 registered), E = 0 vs Meijer G, dps-140 repeat at E = 37, tampers omega*(1 + 1/7) and '
               f'permuted m fail. Drill: omega * 8/7 in spec1_lib (the whole run fails its relation check)')

# drill edits given as (replace, file, None, None) are filled in DRILL_SUBS: first occurrence of `old` -> `new`
DRILL_SUBS = {}
def load_drill_subs():
    p = os.path.join(KIT, 'drills.json')
    if os.path.exists(p):
        DRILL_SUBS.update(json.load(open(p)))

# ---------------------------------------------------------------- integrity
def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''): h.update(b)
    return h.hexdigest()

def _check_sums(p, prefix=''):
    """verify one SHA256SUMS file; paths are relative to its own directory (logs/SHA256SUMS) or to KIT (the others)"""
    bad = []
    base = os.path.dirname(p)
    for line in open(p):
        if not line.strip(): continue
        h, f = line.split(None, 1); f = f.strip()
        q = os.path.join(base if prefix else KIT, f)
        if not os.path.exists(q): bad.append(f'missing {f}')
        elif sha(q) != h: bad.append(f'hash mismatch {f}')
    return bad

def verify_integrity():
    bad = []
    for name in ('SHA256SUMS', 'INPUTS.SHA256'):
        p = os.path.join(KIT, name)
        if not os.path.exists(p): return [f'missing {name}']
        bad += _check_sums(p)
    lp = os.path.join(KIT, 'logs', 'SHA256SUMS')
    if not os.path.exists(lp): bad.append('missing logs/SHA256SUMS')   # v1.2.7 (P1b): the logs manifest is now part of the kit integrity check
    else: bad += _check_sums(lp, prefix='logs/')
    return bad

# ---------------------------------------------------------------- normalisation and criteria
TIMING = [re.compile(r'\(\s*\d+(\.\d+)?\s*s\)'), re.compile(r'\b\d+(\.\d+)?\s?s\)'), re.compile(r'\(\s*\d+(\.\d+)?\s*s\b')]
SKIP = re.compile(r'^(real|user|sys)\s+\d|^(exit|EXIT) \d+$|^\s*$')
def norm(text, mask=()):
    out = []
    for l in text.splitlines():
        if SKIP.match(l): continue
        for r in TIMING: l = r.sub('', l)
        for m in mask: l = re.sub(m, '<numeric diagnostic>', l)
        out.append(l.rstrip())
    return out

def num_diffs(text, tol):
    vals = re.findall(r'diff\s+(-?[0-9.]+e-?\d+|-?[0-9.]+)', text)
    if len(vals) < 5: return False, f'expected >= 5 diff values, found {len(vals)}'
    worst = max(abs(float(v)) for v in vals)
    return worst < tol, f'{len(vals)} diffs, worst {worst:.3g} (tolerance {tol:g})'

# ---------------------------------------------------------------- drill edits
def bump_json(path):
    """add 1/1000 to the first numeric leaf (a number or a rational/expression string) of a JSON file"""
    d = json.load(open(path))
    done = []
    def walk(o):
        if done: return o
        if isinstance(o, dict):
            for k in sorted(o, key=str):
                o[k] = walk(o[k])
                if done: break
            return o
        if isinstance(o, list):
            for i in range(len(o)):
                o[i] = walk(o[i])
                if done: break
            return o
        if isinstance(o, bool): return o
        if isinstance(o, (int, float)):
            done.append(1); return o + 0.001 if isinstance(o, float) else str(Fraction(o) + Fraction(1, 1000))
        if isinstance(o, str) and re.fullmatch(r'\s*-?\d+(/\d+)?\s*', o):
            done.append(1); return str(Fraction(o.strip()) + Fraction(1, 1000))
        if isinstance(o, str) and re.search(r'[0-9]', o) and re.fullmatch(r'[0-9a-zA-Z_+\-*/().^ ]+', o) and len(o) < 20000 and any(c in o for c in '+-*/'):
            done.append(1); return '(' + o + ')+1/1000'
        return o
    d = walk(d)
    if not done: raise RuntimeError(f'bump: no numeric leaf in {path}')
    json.dump(d, open(path, 'w'))

def bump_at(path, keys, rel=None):
    """perturb one chosen data leaf: rationals +1/1000, decimals +1e-20, expressions +1/1000; keys may use '#i' (i-th key/element)"""
    from decimal import Decimal, getcontext
    getcontext().prec = 80
    d = json.load(open(path)); o = d; trail = []
    for k in keys:
        if isinstance(k, str) and k.startswith('#'):
            i = int(k[1:]); k = list(o.keys())[i] if isinstance(o, dict) else i
        trail.append((o, k)); o = o[k]
    parent, k = trail[-1]; v = parent[k]
    if isinstance(v, str) and re.fullmatch(r'\s*-?\d+(/\d+)?\s*', v): nv = str(Fraction(v.strip()) + Fraction(1, 1000))
    elif isinstance(v, str) and re.fullmatch(r'\s*-?\d*\.\d+([eE][-+]?\d+)?\s*', v): nv = str(Decimal(v.strip()) * (1 + Decimal(rel)) if rel else Decimal(v.strip()) + Decimal('1e-20'))
    elif isinstance(v, str): nv = '(' + v + ')+1/1000'
    elif isinstance(v, (int, float)): nv = v + (0.001 if isinstance(v, float) else 0) if isinstance(v, float) else str(Fraction(v) + Fraction(1, 1000))
    else: raise RuntimeError(f'bumpat: unsupported leaf {type(v)} in {path}')
    parent[k] = nv
    json.dump(d, open(path, 'w'))
    return f'{keys} {str(v)[:40]!r} -> {str(nv)[:40]!r}'

def apply_edits(work, edits):
    for e in edits:
        if e[0] == 'bumpat':
            bump_at(os.path.join(work, e[1]), e[2], e[3] if len(e) > 3 else None)
        elif e[0] == 'bump':
            bump_json(os.path.join(work, e[1]))
        elif e[0] == 'rehash':
            f, hf = os.path.join(work, e[1]), os.path.join(work, e[2])
            new = sha(f); base = os.path.basename(e[1]); lines = []
            for l in open(hf):
                if l.strip() and l.split()[-1].endswith(base): l = new + l[64:]
                lines.append(l)
            open(hf, 'w').writelines(lines)
        elif e[0] == 'replace':
            f = os.path.join(work, e[1]); old, new = e[2], e[3]
            if old is None:
                old, new = DRILL_SUBS[e[1]]
            s = open(f).read()
            if old not in s: raise RuntimeError(f'replace: pattern not found in {e[1]}: {old!r}')
            open(f, 'w').write(s.replace(old, new, 1))
        elif e[0] == 'stale':   # v1.2.7 (P1b): set one entry of a SHA256SUMS file to a wrong hash
            hf = os.path.join(work, e[1])
            name = e[2].rsplit('/', 1)[-1]
            lines, done = [], False
            for l in open(hf):
                if not done and l.split()[-1].strip() == name:
                    l = 'd' * 64 + l[64:]; done = True
                lines.append(l)
            if not done: raise RuntimeError(f'stale: entry not found: {name} in {e[1]}')
            open(hf, 'w').writelines(lines)
        elif e[0] == 'delete':   # v1.2.7 (P2): remove one file from a temporary tree copy
            import pathlib
            p = pathlib.Path(os.path.join(work, e[1]))
            if not p.exists(): raise RuntimeError(f'delete: no such file: {e[1]}')
            p.unlink()

# ---------------------------------------------------------------- running
def run_one(name, spec, drill=False, timeout=None, keep=None):
    timeout = timeout or spec.get('timeout', 3600)
    work = tempfile.mkdtemp(prefix='opers-kit-')
    if spec['cwd'] == '.':   # v1.2.7 (P1b): kit-level checks copy the whole kit so manifest drills can tamper with it;
        shutil.rmtree(work)                                  # the entry for this file is fixed up after the copy
        shutil.copytree(KIT, work, symlinks=False, ignore=shutil.ignore_patterns('__pycache__', '.git'))
        _sums = os.path.join(KIT, 'SHA256SUMS')   # v1.2.7: regenerate after removing the self-reference line
        _out = [l for l in open(_sums) if not l.split()[-1].strip().endswith('run_all.py')]
        open(os.path.join(work, 'SHA256SUMS'), 'w').writelines(_out)
    try:
        root = work
        if spec['cwd'] != '.':   # v1.2.7 (P1b): tree checks keep the original layout (fresh copy of tree/ only)
            shutil.copytree(TREE, os.path.join(work, 'tree'), symlinks=False)
            root = os.path.join(work, 'tree')
        argv = spec['argv']
        if drill:
            d = spec['drill']
            if 'edit' in d: apply_edits(root, d['edit'])
            argv = d.get('argv', argv)
        cwd = os.path.join(root, spec['cwd'])
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', **spec['env'])
        t0 = time.time()
        try:
            r = subprocess.run([sys.executable, '-B'] + argv, cwd=cwd, env=env, capture_output=True, text=True, timeout=timeout)
            rc, stdout, stderr = r.returncode, r.stdout, r.stderr
        except subprocess.TimeoutExpired as ex:
            rc, stdout, stderr = 'timeout', (ex.stdout or b'').decode() if isinstance(ex.stdout, bytes) else (ex.stdout or ''), ''
        dt = time.time() - t0
        out = stdout if spec['out'] == 'stdout' else (open(os.path.join(cwd, spec['out'])).read() if os.path.exists(os.path.join(cwd, spec['out'])) else '')
        reasons = []
        if rc != spec['code']: reasons.append(f'exit {rc} (expected {spec["code"]})')
        if spec['ref']:
            ref = open(os.path.join(TREE, spec['cwd'], spec['ref'])).read()
            a, b = norm(out, spec['mask']), norm(ref, spec['mask'])
            if spec['ref_lines'] is not None: b = [b[i] for i in spec['ref_lines']]
            if a != b:
                k = next((i for i in range(min(len(a), len(b))) if a[i] != b[i]), min(len(a), len(b)))
                reasons.append(f'output differs from archived {spec["ref"]} at line {k + 1} ({len(a)} vs {len(b)} lines)')
        if spec['crit']:
            kind, arg = (spec['crit'] + (None,))[:2]
            if kind == 'integrity':   # v1.2.7 (P1b): the gate's own verdict is the criterion
                ok = rc == 0 and 'KIT INTEGRITY FAIL' not in out
                msg = '' if ok else f'integrity gate exited {rc}'
            else:
                ok, msg = num_diffs(out, arg)
            if not ok: reasons.append(msg)
        if spec['post'] and not drill:
            kind, f = spec['post']
            if sha(os.path.join(cwd, f)) != sha(os.path.join(TREE, spec['cwd'], f)): reasons.append(f'{f} differs from the archived file')
        if drill:
            for m in spec['drill'].get('must', []):
                if m not in out + stdout: reasons.append(f'tamper reason {m!r} not in output')
        if keep is not None:
            keep.write(f'\n===== {name}{" [DRILL]" if drill else ""}: argv {argv}  exit {rc}  {dt:.0f}s\n{stdout}')
            if stderr: keep.write(f'--- stderr ---\n{stderr[-4000:]}\n')
        return reasons, dt, out + stdout
    finally:
        shutil.rmtree(work, ignore_errors=True)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true'); ap.add_argument('--full', action='store_true')
    ap.add_argument('--drill', action='store_true'); ap.add_argument('--only'); ap.add_argument('--list', action='store_true')
    ap.add_argument('--log', help='write full outputs to this file')
    a = ap.parse_args()
    load_drill_subs()
    if a.list:
        for n, s in C.items(): print(f'{n:24s} {s["section"]:40s} {"quick" if s["quick"] else "full "}  {s["kind"]}')
        return 0
    bad = verify_integrity()
    if bad:
        print('KIT INTEGRITY FAIL:'); [print('  ' + b) for b in bad[:50]]; return 2
    print('kit integrity: OK (SHA256SUMS, INPUTS.SHA256, logs/SHA256SUMS)')
    names = [n for n, s in C.items() if (a.full or s['quick'])]
    if a.only: names = [n for n in a.only.split(',')]
    keep = open(a.log, 'w') if a.log else None
    fails = 0; T = time.time()
    for n in names:
        s = C[n]
        if a.drill:
            reasons, dt, _ = run_one(n, s, drill=True, keep=keep)
            fired = bool(reasons) and not any(r.startswith('tamper reason') for r in reasons)
            print(f'{"PASS" if fired else "FAIL"}  drill {n:24s} {dt:6.0f}s  {"tamper detected: " + reasons[0] if fired else "TAMPER NOT DETECTED " + "; ".join(reasons)}', flush=True)
        else:
            reasons, dt, _ = run_one(n, s, keep=keep)
            print(f'{"PASS" if not reasons else "FAIL"}  {n:24s} {dt:6.0f}s  [{s["kind"]}; paper: {s["section"]}]{"  " + "; ".join(reasons) if reasons else ""}', flush=True)
            fired = not reasons
        fails += not fired
    print(f'{"ALL PASS" if not fails else f"{fails} FAILED"}  ({len(names)} {"drills" if a.drill else "checks"}, {time.time() - T:.0f}s)')
    return 0 if not fails else 1

if __name__ == '__main__':
    sys.exit(main())
