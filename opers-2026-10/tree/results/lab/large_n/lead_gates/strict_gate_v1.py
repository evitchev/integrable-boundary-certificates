"""Strict gate over the lead's comparison scripts (2026-10-02; Codex review ib-astra-review-14fa9fb finding 4: the comparators printed
mismatches but always exited 0).  The archived comparators are left untouched.  For each one this gate re-runs it on its archived prediction
(or on an override given as NAME=path), from the repo root resolved from __file__, and FAILS (exit 2) if
  - a comparison line reports a non-zero mismatch count,
  - a negative-control line reports zero mismatches,
  - a line is neither, unless it is one of the enumerated exceptional lines below (known, recorded cells),
  - the run exits non-zero or prints nothing,
  - (archived prediction only) the output differs from the archived log.
Usage: python3 strict_gate.py [NAME=prediction.json ...]"""
import re, subprocess, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[4]
L = 'results/lab/large_n/'
JOBS = {'VIR7': ('lead_vir7_check/compare_vir7.py', 'lead_vir7_check/PREDICTION_VIR7.json', 'lead_vir7_check/compare_vir7.log'),
        'VIR4D': ('lead_vir4d_check/compare_vir4d.py', 'vir4d_prereg/PREDICTION_VIR4d.json', 'lead_vir4d_check/compare_vir4d.log'),
        'VIR5A': ('lead_vir5a_check/compare_vir5a.py', 'lead_vir5a_check/PREDICTION_VIR5a.json', 'lead_vir5a_check/compare_vir5a.log'),
        'VIR5C': ('lead_vir5c_check/compare_spin15.py', 'lead_vir5c_check/PREDICTION_VIR5c.json', 'lead_vir5c_check/compare_spin15.log'),
        'SOL2F': ('lead_sol2f_check/compare_sol2f.py', 'lead_sol2f_check/PREDICTION_SOL2F.json', 'lead_sol2f_check/compare_sol2f.log')}
EXCEPTIONAL = {  # recorded rank-drop cells: the oper's spin-11 coefficient vanishes there by the Beta-pole mechanism (items 303/307)
    'k=4,t=7/3 M=1/4 predicted spin 11 = 0; certified table at t: NONZERO (28 terms)',
    'k=4,t=7/3 M=-1/5 predicted spin 11 = 0; certified table at t: NONZERO (28 terms)',
    '--- negative controls ---'}
TIMING = re.compile(r'^(real|user|sys)\s')
def classify(line):
    if line in EXCEPTIONAL: return None
    neg = line.startswith('NEG')
    if neg:
        m = re.search(r'mismatches (\d+)', line); return None if m and int(m.group(1)) > 0 else 'negative control did not fire'
    errs = []
    for m in re.finditer(r'(\d+):(\d+)/\d+\(neg (\d+)\)', line):
        if int(m.group(2)): errs.append(f'spin {m.group(1)} mismatch {m.group(2)}')
        if not int(m.group(3)): errs.append(f'spin {m.group(1)} negative control did not fire')
    for m in re.finditer(r'(\d+):(\d+)/\d+(?!\(|\d)', line):
        if int(m.group(2)): errs.append(f'spin {m.group(1)} mismatch {m.group(2)}')
    m = re.search(r'mismatches (\d+)', line)
    if m and int(m.group(1)): errs.append(f'mismatch {m.group(1)}')
    m = re.search(r'NEG[^:]*: (\d+)', line)
    if m and not int(m.group(1)): errs.append('negative control did not fire')
    recognised = bool(re.search(r'mismatches \d+|\d+:\d+/\d+|:no-charge', line))
    if not recognised: errs.append('unrecognised line')
    return '; '.join(errs) or None
over = dict(a.split('=', 1) for a in sys.argv[1:])
FAIL = []
for name, (script, pred, log) in JOBS.items():
    p = over.get(name, L + pred)
    r = subprocess.run([sys.executable, '-B', L + script, p], cwd=ROOT, capture_output=True, text=True)
    out = [x for x in r.stdout.splitlines() if x.strip()]
    if r.returncode or not out: FAIL.append(f'{name}: exit {r.returncode}, {len(out)} lines, stderr {r.stderr[-300:]!r}')
    for x in out:
        e = classify(x)
        if e: FAIL.append(f'{name}: {e} :: {x}')
    if name not in over:
        arch = [x for x in (ROOT / L / log).read_text().splitlines() if x.strip() and not TIMING.match(x)]
        if arch != out: FAIL.append(f'{name}: output differs from archived log')
    print(f'{name}: {len(out)} lines, exit {r.returncode}', flush=True)
for f in FAIL: print('FAIL', f)
print('GATE', 'FAIL' if FAIL else 'PASS')
sys.exit(2 if FAIL else 0)
