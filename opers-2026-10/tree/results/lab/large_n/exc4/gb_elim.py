"""EXC4 B: from msolve's grevlex Groebner basis (h_fourq_<TAG>_four.msgb) to the lex basis by FGLM in Singular (zero-dimensional), then the
univariate eliminants in ee (monic I_3 over the sigma-invariant four-point solutions) and ii (monic I_1), and vdim.  Usage: gb_elim.py <TAG>"""
import re, subprocess, sys, time
sys.dont_write_bytecode = True
TAG = sys.argv[1]
s = open('h_fourq_%s_four.msgb' % TAG).read()
vars_ = re.search(r'#variable order:\s*(.*)', s).group(1).replace(' ', '')
body = s[s.index('[') + 1:s.rindex(']')]
polys = [p.strip() for p in re.split(r',\s*\n', body) if p.strip()]
t0 = time.time()
script = ['ring R = 0,(%s),dp;' % vars_, 'ideal G =', ',\n'.join(polys) + ';', 'attrib(G,"isSB",1);', 'print("vdim:"); vdim(G);',
          'ring Rl = 0,(%s),lp;' % vars_, 'ideal Gl = fglm(R, G);', 'print("lexsize:"); size(Gl);',
          'int i; ideal V; for (i = 1; i <= size(Gl); i++) { V = variables(Gl[i]); if (size(V) == 1) { if (V[1] == ee) { print("eliminant ee:"); print(Gl[i]); } if (V[1] == ii) { print("eliminant ii:"); print(Gl[i]); } if (V[1] == z) { print("eliminant z:"); print(Gl[i]); } } }',
          'quit;']
fn = 'gb_elim_%s.sing' % TAG
open(fn, 'w').write('\n'.join(script) + '\n')
res = subprocess.run(['<home>/ib-lab/solvers/env/bin/Singular', '-q'], stdin=open(fn), capture_output=True, text=True, timeout=7000).stdout
open(fn.replace('.sing', '.out'), 'w').write(res)
print(TAG, '(%d GB elements, %.0fs)' % (len(polys), time.time() - t0)); print(res[:1500])
