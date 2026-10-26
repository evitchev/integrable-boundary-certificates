#!/usr/bin/env python3
"""v1.2.7 (P1b): kit integrity gate as a standalone check target.

Runs the same verify_integrity() that opers-2026-10/run_all.py runs before any
check -- now including opers-2026-10/logs/SHA256SUMS -- and exits non-zero with
a 'KIT INTEGRITY FAIL' report if anything is missing or stale.
"""
import os, sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'opers-2026-10'))
import run_all

bad = run_all.verify_integrity()
if bad:
    print('KIT INTEGRITY FAIL:')
    for b in bad[:50]:
        print('  ' + b)
    sys.exit(1)
print('kit integrity: OK (SHA256SUMS, INPUTS.SHA256, logs/SHA256SUMS)')
