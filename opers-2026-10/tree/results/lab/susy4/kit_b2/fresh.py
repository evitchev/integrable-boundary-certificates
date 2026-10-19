"""kit_b2/fresh.py -- kit helper (not an archived file) for the B2/B2b checks of kit v1.2.6.

Usage (cwd = the stage directory):
    python3 -B ../../kit_b2/fresh.py [--rm FILE ...] [--same-json FILE ...] [--ignore-key KEY ...] -- SCRIPT [ARGS ...]
--rm FILE         delete FILE (an output of SCRIPT that is also archived in the stage directory) before the run.
--same-json FILE  read the archived FILE, delete it, run SCRIPT, then require that SCRIPT rewrote FILE and that it equals the
                  archived one as a JSON object (key order not significant; keys named by --ignore-key, e.g. timings, dropped
                  at every level on both sides).  On a mismatch the helper exits with status 3 (the check expects 0).
SCRIPT is the archived script, run unchanged as __main__ with sys.argv = [SCRIPT, ARGS ...] and its own directory first on
sys.path, exactly as `python3 SCRIPT ARGS` would.  Purpose: a check that compares SCRIPT's written output with the archived
file must not pass because the archived copy was left in place.  The helper writes only to stderr, so SCRIPT's stdout
is compared with the archived log unchanged.
"""
import json
import os
import runpy
import sys

sys.dont_write_bytecode = True
args = sys.argv[1:]
if '--' not in args:
    sys.exit('usage: fresh.py [--rm FILE ...] [--same-json FILE ...] [--ignore-key KEY ...] -- SCRIPT [ARGS ...]')
k = args.index('--')
opts, cmd = args[:k], args[k + 1:]


def opt(name):
    return [opts[i + 1] for i in range(len(opts)) if opts[i] == name]


ignore = set(opt('--ignore-key'))


def strip(o):
    if isinstance(o, dict):
        return {a: strip(b) for a, b in o.items() if a not in ignore}
    if isinstance(o, list):
        return [strip(x) for x in o]
    return o


keep = {}
for f in opt('--same-json'):
    keep[f] = strip(json.load(open(f)))
for f in opt('--rm') + opt('--same-json'):
    if os.path.exists(f):
        os.remove(f)
        print(f'[kit_b2/fresh] removed archived copy of {f} before the run', file=sys.stderr, flush=True)
script = cmd[0]
sys.argv = cmd
sys.path.insert(0, os.path.dirname(os.path.abspath(script)))
rc = 0
try:
    runpy.run_path(script, run_name='__main__')
except SystemExit as ex:          # the script's own sys.exit(n): compare outputs first, then exit with its status
    rc = ex.code if isinstance(ex.code, int) else (0 if ex.code is None else 1)
sys.stdout.flush()
bad = []
for f, ref in keep.items():
    if not os.path.exists(f):
        bad.append(f'{f} not rewritten')
    elif strip(json.load(open(f))) != ref:
        bad.append(f'{f} differs from the archived file as a JSON object')
    else:
        print(f'[kit_b2/fresh] {f} equals the archived file as a JSON object' + (f' (ignoring {sorted(ignore)})' if ignore else ''),
              file=sys.stderr, flush=True)
if bad:
    print('[kit_b2/fresh] FAIL: ' + '; '.join(bad), file=sys.stderr, flush=True)
    sys.exit(3)
sys.exit(rc)
