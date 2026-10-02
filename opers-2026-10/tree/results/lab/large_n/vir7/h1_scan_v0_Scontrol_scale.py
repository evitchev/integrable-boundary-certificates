"""VIR7 T0/T1: scan of class U (Sol 1) at first quantum order with formula G1, plus controls.
Exit 0 = a member passes T1 and controls behave; 3 = no member passes, controls behave; 1 = a control misbehaves.
Mode: real | tamper (Codex A_(1,1) + 1/1000 on the Sol-2 plant, after the seal-hash guard)."""
import hashlib, json, sys
import sympy as sp
sys.dont_write_bytecode = True
import g1_lib as G
SEAL = open('SEAL_VIR7.sha256').read().split()[0]
assert hashlib.sha256(open('SEAL_VIR7.md', 'rb').read()).hexdigest() == SEAL
mode = sys.argv[1] if len(sys.argv) > 1 else 'real'
KMAX = 20
AUD = '<repo>/results/lab/audits/codex_wardsol12b_27403fc_2026-09-28/'
ks, ts = sp.symbols('k t')
l0, l1, ZERO = G.l0, G.l1, G.ZERO


def load_law(sol):
    raw = open(AUD + 'AMPLITUDES_sol%d_loss1.json' % sol, 'rb').read()
    law = json.loads(raw)
    f = lambda d: sp.sympify(d['numerator'], locals={'k': ks, 't': ts}, rational=True) / sp.sympify(d['denominator'], locals={'k': ks, 't': ts}, rational=True)
    return {int(j): f(v) for j, v in law['shift_amplitudes'].items()}, f(law['p']), f(law['q']), f(law['r']), hashlib.sha256(raw).hexdigest()


def Kd(d, p, q, r, a):
    return sp.binomial(d, a) * sp.rf(p, a) * sp.rf(q, d - a) * r ** (d - a) / sp.rf(p, d)


def ward(law, K, t, tamper=False):
    amps, p, q, r, _ = law
    sub = {ks: K, ts: t}
    p, q, r = [sp.Rational(x.subs(sub)) for x in (p, q, r)]
    vals = {j: sp.Rational(v.subs(sub)) for j, v in amps.items()}
    if tamper:
        vals[1] = vals[1] + sp.Rational(1, 1000)
    lay = {(a, K - 1 - a): sum(v * Kd(K - 1, p + j, q, r, a) for j, v in vals.items()) for a in range(K)}
    top = {(a, K - a): Kd(K, p, q, r, a) for a in range(K + 1)}
    return lay, top


def member(name, t):
    """returns dict(n, M, factors, strings, alpha, beta)"""
    k = 2 * (t - 1) / (3 - t)
    def U(beta_P, beta_L, a):
        n = 2 * beta_P + a * (2 + 2 * beta_L - 2)
        fac = [(l1, beta_P + a * beta_L), (-l1, beta_P + a * beta_L), (l0, a), (-l0, a), (ZERO, -2 * a)]
        strs = [(l0, a, 1), (-l0, a, 1), (ZERO, a, -2)] + ([(l1, a, beta_L), (-l1, a, beta_L)] if beta_L else [])
        return dict(n=n, M=1 / a, factors=fac, strings=strs, alpha=a, beta=beta_P + a * beta_L)
    if name.startswith('U1_'):
        u = int(name[3:]); return U(1, 3 - u, k / (u * k + 4))
    if name.startswith('U2_'):
        u = int(name[3:]); return U(2, 3 - 2 * u, k / (u * k + 2))
    if name == 'U4':
        return U(4, 3, k)
    if name == 'Gcc':          # cc's beta = 1 symbol: one frozen block of length -2 alpha
        a = (t - 1) / (t + 3)
        return dict(n=sp.Integer(2), M=1 / a, factors=[(l1, 1), (-l1, 1), (l0, a), (-l0, a), (ZERO, -2 * a)],
                    strings=[(l0, a, 1), (-l0, a, 1), (ZERO, -2 * a, 1)], alpha=a, beta=sp.Integer(1))
    if name == 'S':            # single-string form of Suzuki + centrifugal
        a = (t - 1) / 2
        return dict(n=2 + a, M=1 / a, factors=[(l1, 1), (-l1, 1), (l0, a)], strings=[(l0, a, 1)], alpha=None, beta=sp.Integer(1))
    if name in ('SOL2', 'SOL2_tamper'):
        a = k + (sp.Rational(1, 10) if name == 'SOL2_tamper' else 0)
        return dict(n=2 * k + 3, M=1 / k, factors=[(l0, k), (-l0, k), (l1, 1), (-l1, 1), (ZERO, 1)], strings=[(l0, a, 1), (-l0, a, 1)], alpha=k, beta=sp.Integer(1))
    raise ValueError(name)


def strings_expand(strs):
    out = []
    for e, a, mult in strs:
        # mult copies of a string of length a (mult < 0: inverse strings, c_1 enters with the opposite sign)
        out.append((e, a, mult))
    return out


def run_member(name, law, fibres, tamper=False, verbose=True):
    """returns (top_ok, t1_ok, info)"""
    top_ok = t1_ok = True
    info = []
    for tstr in fibres:
        t = sp.Rational(tstr)
        rho2 = 2 / (t * t - 1)
        m = member(name, t)
        n, M = sp.Rational(m['n']), sp.Rational(m['M'])
        c2 = n * n * (M + 1)
        # b(z): strings with multiplicity
        N = 2 * KMAX
        sigma2 = (n * M) ** 2
        extra = [ZERO] * (N + 1)
        for e, a, mult in m['strings']:
            a = sp.Rational(a)
            c1 = -a * (a * a - 1) / 24
            g2 = G.geom(e, N, 2)
            extra = [extra[i] + g2[i] * G.q_(c1 * sigma2 * mult) for i in range(N + 1)]
        tb = tt = fb = ft = 0
        invC = None
        witness = None
        for K in range(2, KMAX + 1):
            try:
                lay, top = ward(law, K, t, tamper)
            except ZeroDivisionError:
                continue
            if any(v.has(sp.zoo, sp.nan) for v in list(lay.values()) + list(top.values())):
                continue
            topP, l1P = G.charge(2 * K, n, M, m['factors'], [], extra_b=extra[:2 * K + 1])
            if name == 'S':
                # X scale from the top at K = 2 (control only): try the two signs of CX = +-CY
                CY = sp.Integer(1); CX = None
                for cand in (sp.Integer(1), sp.Integer(-1), sp.Rational(1, 2), sp.Integer(2), sp.Rational(1, 4), sp.Integer(4)):
                    r0 = G.to_record(G.charge(4, n, M, m['factors'], [])[0], ZERO, 2, cand, CY, 0)
                    w2 = ward(law, 2, t)[1]
                    if r0 is not None and all(r0[1].get(key, 0) == v for key, v in w2.items()):
                        CX = cand; break
                if CX is None:
                    return False, False, ['S: no X scale in the tried set reproduces the K = 2 top at t = %s' % tstr]
            else:
                CX, CY = 1 / sp.Rational(m['alpha']), 1 / sp.Rational(m['beta'])
            # C-independent part (rho^2 shift of the top) and the 1/C part
            r_shift = G.to_record(topP, ZERO, K, CX, CY, rho2)
            r_full = G.to_record(topP, l1P, K, CX, CY, rho2)
            if r_shift is None:
                continue
            for key, v in top.items():
                tt += 1; tb += r_shift[1].get(key, 0) != v
            base = r_shift[0]; one = {key: r_full[0].get(key, 0) - base.get(key, 0) for key in set(r_full[0]) | set(base)}
            # record layer = base + (1/C) * one   (momenta l^2 = C * (CX PX^2, CY pi^2))
            for key in sorted(lay):
                tgt = lay[key] - base.get(key, 0)
                if invC is None:
                    if one.get(key, 0) != 0:
                        invC = tgt / one[key]
                    elif tgt != 0:
                        ft += 1; fb += 1
                    continue
                ft += 1
                if invC * one.get(key, 0) != tgt:
                    fb += 1
                    if witness is None:
                        witness = (K, key)
        top_ok &= tb == 0
        t1_ok &= (fb == 0 and ft > 0)
        rule = None if invC in (None, 0) else (1 / invC == c2)
        info.append('t = %-4s n = %-7s M = %-7s top mismatches %d/%d; loss-1 mismatches %d/%d; fitted C = %s; C == n^2(M+1) = %s: %s; first witness %s'
                    % (tstr, n, M, tb, tt, fb, ft, None if invC in (None, 0) else 1 / invC, c2, rule, witness))
    if verbose:
        print('%-11s T0 %s  T1 %s' % (name, 'PASS' if top_ok else 'FAIL', 'PASS' if t1_ok else 'FAIL'))
        for line in info:
            print('     ' + line)
        sys.stdout.flush()
    return top_ok, t1_ok, info


FIB = ['2', '9/4', '3/2', '4', '-2']
law1 = load_law(1)
law2 = load_law(2)
print('Sol-1 loss-1 input sha256', law1[4]); print('Sol-2 loss-1 input sha256', law2[4])
status_ctrl = True
if mode == 'tamper':
    _, ok, _ = run_member('SOL2', law2, FIB, tamper=True)
    print('scoring tamper on the Sol-2 plant:', 'fires' if not ok else 'DOES NOT FIRE')
    sys.exit(0 if not ok else 1)
print('--- controls')
_, ok, _ = run_member('SOL2', law2, FIB); status_ctrl &= ok
_, ok, _ = run_member('SOL2_tamper', law2, FIB); status_ctrl &= (not ok)
tS, okS, _ = run_member('S', law1, FIB)
tG, okG, _ = run_member('Gcc', law1, FIB)
print('--- class U')
passed = []
for name in ('U1_3', 'U1_2', 'U1_1', 'U1_0', 'U2_1', 'U2_0', 'U4'):
    t0, ok, _ = run_member(name, law1, FIB)
    if t0 and ok:
        passed.append(name)
print('members passing T0 and T1:', passed)
print('controls: Sol-2 plant/tamper behave: %s; S fails T1: %s (top %s); Gcc fails T1: %s (top %s)' % (status_ctrl, not okS, tS, not okG, tG))
sys.exit(1 if not status_ctrl else (0 if passed else 3))
