#!/usr/bin/env python3
"""Re-run the N = 1 super-KdV cylinder eigenvalue computation and compare with the archived results.

Usage:  python3 run_checks.py [--quick] [--drill] [--log FILE]

  (default)  check 1 (evaluator validation, ~15 s) and check 2 (tower and eigenvalues, ~40 min)
  --quick    check 1 only
  --drill    the same checks with a deliberate tamper (the Ramond fermion kernel 1/(e^eps - 1) - 1/eps gets an extra
             eps/100); each tampered run must FAIL

Every path is resolved from this file's location.  The files are first checked against SHA256SUMS.  Each check runs
in a fresh temporary copy of code/ and results/, so the archived files are never overwritten.
  check 1: code/super_cyl_vev_validate.py; its output, with run times removed, must equal logs/super_cyl_vev_validate.log.
  check 2: code/super_cyl_tower_fit.py 20; the file it writes, results/lab/susy4/tanabe/computed_eigenvalues.json,
           must be byte-identical to the archived one (and the .txt likewise).
Exit status: 0 iff every requested check (or drill) passes; 1 otherwise; 2 on an integrity error.
"""
import argparse, hashlib, os, re, shutil, subprocess, sys, tempfile, time

HERE = os.path.dirname(os.path.abspath(__file__))
TAMPER = ('code/super_cyl_vev.py', '1 / (sp.exp(_eps) - 1) - 1 / _eps, n)', '1 / (sp.exp(_eps) - 1) - 1 / _eps + _eps / 100, n)')
OUTS = ('results/lab/susy4/tanabe/computed_eigenvalues.json', 'results/lab/susy4/tanabe/computed_eigenvalues.txt')


def sha(p):
    with open(p, 'rb') as f:
        return hashlib.sha256(f.read()).hexdigest()


def integrity():
    bad = []
    for line in open(os.path.join(HERE, 'SHA256SUMS')):
        if not line.strip():
            continue
        h, f = line.split(None, 1)
        p = os.path.join(HERE, f.strip())
        if not os.path.exists(p):
            bad.append('missing ' + f.strip())
        elif sha(p) != h:
            bad.append('hash mismatch ' + f.strip())
    return bad


TIMING = re.compile(r'\(\s*\d+(\.\d+)?\s*s\)|\b\d+\.\d+\s*s\b')


def norm(text):
    return [TIMING.sub('', l).rstrip() for l in text.splitlines() if l.strip()]


def run(argv, drill, timeout):
    work = tempfile.mkdtemp(prefix='super-kdv-')
    try:
        for d in ('code', 'results'):
            shutil.copytree(os.path.join(HERE, d), os.path.join(work, d))
        if drill:
            f, old, new = TAMPER
            p = os.path.join(work, f)
            s = open(p).read()
            if old not in s:
                raise RuntimeError('tamper pattern not found')
            open(p, 'w').write(s.replace(old, new, 1))
        for f in OUTS:
            os.remove(os.path.join(work, f))
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
        t0 = time.time()
        try:
            r = subprocess.run([sys.executable, '-B'] + argv, cwd=os.path.join(work, 'code'), env=env,
                               capture_output=True, text=True, timeout=timeout)
            rc, out = r.returncode, r.stdout + r.stderr
        except subprocess.TimeoutExpired:
            rc, out = 'timeout', ''
        written = {f: (sha(os.path.join(work, f)) if os.path.exists(os.path.join(work, f)) else None) for f in OUTS}
        return rc, out, written, time.time() - t0
    finally:
        shutil.rmtree(work, ignore_errors=True)


def check1(drill):
    rc, out, _, dt = run(['super_cyl_vev_validate.py'], drill, 1800)
    ref = open(os.path.join(HERE, 'logs', 'super_cyl_vev_validate.log')).read()
    reasons = []
    if rc != 0:
        reasons.append(f'exit {rc}')
    if norm(out) != norm(ref):
        reasons.append('output differs from logs/super_cyl_vev_validate.log')
    return reasons, dt, out


def check2(drill):
    rc, out, written, dt = run(['super_cyl_tower_fit.py', '20'], drill, 3 * 3600)
    reasons = []
    if rc != 0:
        reasons.append(f'exit {rc}')
    for f, h in written.items():
        if h != sha(os.path.join(HERE, f)):
            reasons.append(f'{os.path.basename(f)} {"not written" if h is None else "differs from the archived file"}')
    return reasons, dt, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--drill', action='store_true')
    ap.add_argument('--log')
    a = ap.parse_args()
    bad = integrity()
    if bad:
        print('INTEGRITY FAIL:')
        for b in bad:
            print('  ' + b)
        return 2
    print('integrity: OK (SHA256SUMS)', flush=True)
    keep = open(a.log, 'w') if a.log else None
    checks = [('evaluator_validation', check1)] + ([] if a.quick else [('tower_eigenvalues', check2)])
    fails = 0
    for name, fn in checks:
        reasons, dt, out = fn(a.drill)
        if keep:
            keep.write(f'\n===== {name}{" [DRILL]" if a.drill else ""}  {dt:.0f}s\n{out}')
        if a.drill:
            ok = bool(reasons)
            print(f'{"PASS" if ok else "FAIL"}  drill {name:22s} {dt:6.0f}s  '
                  f'{"tamper detected: " + reasons[0] if ok else "TAMPER NOT DETECTED"}', flush=True)
        else:
            ok = not reasons
            print(f'{"PASS" if ok else "FAIL"}  {name:22s} {dt:6.0f}s  {"; ".join(reasons)}', flush=True)
        fails += not ok
    print('ALL PASS' if not fails else f'{fails} FAILED')
    return 0 if not fails else 1


if __name__ == '__main__':
    sys.exit(main())
