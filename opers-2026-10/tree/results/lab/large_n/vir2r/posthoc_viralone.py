"""VIR2-R POST-HOC (informational, labelled): commutant dimension mod d with ONLY the Virasoro screening exp(p phi) (no X-screening pair), and with
exp(p phi) + ONE member V_+ only, at Sol 1 p = 8/5, weights 4, 6, 8 -- shows the pair is what cuts the commutant to 1."""
import sys, os
sys.argv = ['x', '--nocomm']
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
src = open(os.path.join(HERE, 'vir2r.py')).read()
pre = src[:src.index('# ======================= run')]
pre = pre.replace("LOG = open(os.path.join(HERE, 'run_VIR2R.log'), 'w')", "LOG = open(os.path.join(HERE, 'run_posthoc_viralone.inner.log'), 'w')")
cm = src[src.index('def commutant_dim('):src.index('if not NOCOMM:')]
g = {'__file__': os.path.join(HERE, 'vir2r.py'), '__name__': 'ph'}; exec(pre + cm, g)
F, QS, XC, FC = g['F'], g['QS'], g['XC'], g['FC']
p = F(8, 5); a0 = F(3, 5); beta = -p / 2; QS.d = F(1)
for W in (4, 6, 8):
    d1, _, _ = g['commutant_dim'](W, [{XC: QS(0), FC: QS(p)}])
    d2, _, _ = g['commutant_dim'](W, [{XC: QS(0), FC: QS(p)}, {XC: QS(a0), FC: QS(beta)}])
    print(f'Sol 1 p = 8/5, weight {W}: commutant mod d with exp(p phi) only: {d1}; with exp(p phi) + V_+ only: {d2}')
