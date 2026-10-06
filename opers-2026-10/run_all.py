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
def check(name, section, kind, cwd, argv, ref, out='stdout', code=0, quick=True, env=None, crit=None, post=None, drill=None, note='', mask=None, ref_lines=None):
    C[name] = dict(section=section, kind=kind, cwd=cwd, argv=argv, ref=ref, out=out, code=code, quick=quick,
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

def verify_integrity():
    bad = []
    for name in ('SHA256SUMS', 'INPUTS.SHA256'):
        p = os.path.join(KIT, name)
        if not os.path.exists(p): return [f'missing {name}']
        for line in open(p):
            if not line.strip(): continue
            h, f = line.split(None, 1); f = f.strip()
            q = os.path.join(KIT, f)
            if not os.path.exists(q): bad.append(f'missing {f}')
            elif sha(q) != h: bad.append(f'hash mismatch {f}')
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

# ---------------------------------------------------------------- running
def run_one(name, spec, drill=False, timeout=3600, keep=None):
    work = tempfile.mkdtemp(prefix='opers-kit-')
    try:
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
            ok, msg = num_diffs(out, spec['crit'][1])
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
    print('kit integrity: OK (SHA256SUMS, INPUTS.SHA256)')
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
