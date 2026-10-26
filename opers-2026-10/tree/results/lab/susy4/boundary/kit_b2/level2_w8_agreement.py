#!/usr/bin/env python3
"""v1.2.7 (P2): agreement of the B2 w8 real-projection row with the rational membership test.

The level2.py CONTROLS block prints, at every fibre:

    real projection at s in the (X, phi) frame: D_0 in span(I_3^cert) + Im d: <value>  (must be True)

The label "(must be True)" is a pre-result expectation; the row is the real-dictionary
projection at w2 = 8, printed for consistency with the rational membership test (the
"(q^k, weight-w) pieces in Solution 1's commutant ... ALL: ..." line of the same log).
This check enforces exactly that agreement on the archived run: the row must read
False and the same log's w2 = 8 membership line must end "ALL: False".
Usage: level2_w8_agreement.py Z|K
"""
import os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
LOGS = os.path.join(HERE, '..', 'b2', 'logs')
ROW = re.compile(r'real projection at s in the \(X, phi\) frame: D_0 in span\(I_3\^cert\) \+ Im d: (True|False)\s*\(must be True\)')

def main():
    tower = sys.argv[1] if len(sys.argv) > 1 else 'Z'
    fails = []
    logs = sorted(os.path.join(LOGS, f) for f in os.listdir(LOGS)
                  if re.fullmatch(r'level2_' + tower + r'_t\d+o?\d*\.log', f))
    if not logs:
        print(f'AGREEMENT CHECK FAILED: no level2_{tower}_*.log under {LOGS}')
        return 1
    for lg in logs:
        text = open(lg).read()
        lines = text.splitlines()
        mrow = [i for i, l in enumerate(lines) if ROW.search(l)]
        if len(mrow) != 1:
            fails.append(f'{os.path.basename(lg)}: expected exactly one w8 row, found {len(mrow)}')
            continue
        val = ROW.search(lines[mrow[0]]).group(1)
        # the same log's w2 = 8 membership line (first "pieces in Solution 1's commutant" after "== w2 = 8")
        sec = None
        for i, l in enumerate(lines):
            if l.startswith('== w2 = 8'):
                for j in range(i + 1, len(lines)):
                    if "pieces in Solution 1's commutant" in lines[j]:
                        sec = lines[j]; break
                break
        if sec is None:
            fails.append(f'{os.path.basename(lg)}: no w2 = 8 membership line')
            continue
        tail = 'ALL: False' if sec.rstrip().endswith('ALL: False') else ('ALL: True' if sec.rstrip().endswith('ALL: True') else 'UNPARSEABLE')
        print(f'  {os.path.basename(lg)}: w8 row = {val}; w2 = 8 membership = {tail}')
        if val != 'False' or not sec.rstrip().endswith('ALL: False'):
            fails.append(f'{os.path.basename(lg)}: disagreement (row {val}, membership {tail})')
    if fails:
        print('AGREEMENT CHECK FAILED:')
        [print('  ' + f) for f in fails]
        return 1
    print(f'w8 real-projection row agrees with the rational membership test in {len(logs)} level2_{tower} logs')
    return 0

if __name__ == '__main__':
    sys.exit(main())
