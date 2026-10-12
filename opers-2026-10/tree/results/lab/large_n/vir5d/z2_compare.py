"""VIR5d Z2: compare the certified kernel at the drop points with the deformation derivatives of the operator.
Vectors = coefficient vectors over the monomials PX^(2a) pi^(2b), a + b <= K.  All ranks exact (sympy Rational).
Modes: real | tamper (one coefficient of a kernel class + 1/1000, after the hash guard).
Exit (real): 0 = H-A passes at every 2-dim-kernel point and controls fire; 2 = H-A fails as sealed; 1 = control failure."""
import hashlib, json, pickle, sys
import sympy as sp
sys.dont_write_bytecode = True
SEAL = open('SEAL_VIR5d.sha256').read().split()[0]
assert hashlib.sha256(open('SEAL_VIR5d.md', 'rb').read()).hexdigest() == SEAL
mode = sys.argv[1] if len(sys.argv) > 1 else 'real'
ROOT = '<repo>/results/lab/'
t, p, PX, PI = sp.symbols('t p PX PI')
Fsym = {k: sp.Poly(sp.sympify(v), PX, PI) for k, v in pickle.load(open('cft_F.pkl', 'rb')).items()}


def table_poly(sol, s, tv):
    """certified vacuum eigenvalue from the record's tables, in (PX, pi), un-normalised"""
    rho2 = sp.Rational(2) / (tv**2 - 1)
    tab = json.load(open(ROOT + 'vev/vev_sol%d_w%d.json' % (sol, s + 1)))
    expr = 0
    for key, val in tab['coefficients'].items():
        a_, b_ = [int(x.split('^')[1]) for x in key.split()]
        f = sp.sympify(val, locals={'t': t}, rational=True)
        assert sp.fraction(sp.together(f))[1].subs(t, tv) != 0
        expr += f.subs(t, tv) * (-PX**2) ** (a_ // 2) * (rho2 - PI**2) ** (b_ // 2)
    return sp.Poly(sp.expand(expr), PX, PI)


def vec(d, mons):
    return [sp.Rational(d.get(m, 0)) for m in mons]


def to_dict(P):
    return {(a // 2, b // 2): c for (a, b), c in P.terms()}


def load_keys(d):
    return {(int(k.split(',')[0]) // 2, int(k.split(',')[1]) // 2): sp.Rational(v) for k, v in d.items()}


def swap(d):
    return {(b, a): c for (a, b), c in d.items()}


def rank(rows):
    return sp.Matrix(rows).rank() if rows else 0


def show(name, d, K):
    out = []
    for l_ in range(0, K + 1):
        row = ['%s X^%d Y^%d' % (d[(a, K - l_ - a)], a, K - l_ - a) for a in range(K - l_, -1, -1) if d.get((a, K - l_ - a), 0) != 0]
        if row:
            out.append('      loss %d: ' % l_ + ' ; '.join(row))
    print('   %s  (X = PX^2, Y = pi^2):' % name)
    print('\n'.join(out) if out else '      0')


status = 0
POINTS = [('5/3', '1', 3), ('7/5', '1/2', 5), ('11/5', '3', 4), ('17/7', '5', 5)]
for tstr, kstr, K in POINTS:
    tv = sp.Rational(tstr)
    mons = [(a, b) for a in range(K + 1) for b in range(K + 1 - a)]
    kd = json.load(open('z2_kernel_t%s_K%d.json' % (tstr.replace('/', '_'), K)))
    ker = [load_keys(c['eigenvalue']) for c in kd['classes'] if c['in_HV']]
    assert len(ker) == kd['kernel_dimension']
    if mode == 'tamper':
        ker[-1] = dict(ker[-1]); ker[-1][(1, 1)] = ker[-1].get((1, 1), 0) + sp.Rational(1, 1000)
    print('\n== t = %s (k = %s), K = %d (spin %d): certified kernel dimension %d' % (tstr, kstr, K, 2 * K - 1, len(ker)))
    KR = [vec(e, mons) for e in ker]
    rk = rank(KR)
    assert rk == len(ker)
    # gate: the certified Sol-3 table charge lies in the kernel
    S = to_dict(table_poly(3, 2 * K - 1, tv))
    inS = rank(KR + [vec(S, mons)]) == rk
    print('   gate: Sol-3 table charge at this fibre lies in the computed kernel: %s' % inS)
    if not inS and mode == 'real':
        status = 1
    if tstr == '5/3':
        S2 = to_dict(table_poly(2, 5, tv))
        in2 = rank(KR + [vec(S2, mons)]) == rk
        ind = rank([vec(S, mons), vec(S2, mons)]) == 2
        print('   handle: Sol-2 table charge at t = 5/3 lies in the kernel: %s; independent of the Sol-3 charge: %s' % (in2, ind))
    # (ii) swap invariance
    sw = all(rank(KR + [vec(swap(e), mons)]) == rk for e in ker)
    print('   (ii) kernel invariant under PX <-> pi: %s' % sw)
    Ssym = swap(S) == S
    print('        Sol-3 charge symmetric: %s' % Ssym)
    # operator side
    for da, db, strp, name in (('1', 'm1', False, 'A  (antisymmetric power deformation)'), ('1', '1', False, 'S~ (symmetric power deformation, fixed n)'),
                               ('1', '0', False, 'A_X (X-power only)'), ('1', 'm1', True, 'A stripped')):
        fn = 'z2_def_k%s_K%d_%s_%s%s.json' % (kstr.replace('/', '_'), K, da, db, '_stripped' if strp else '')
        try:
            dd = json.load(open(fn))
        except FileNotFoundError:
            continue
        R1 = load_keys(dd['R1'])
        zero0 = dd['R0_identically_zero']
        inK = rank(KR + [vec(R1, mons)]) == rk if zero0 else None
        anti = all(R1.get((b, a), 0) == -c for (a, b), c in R1.items())
        symm = swap(R1) == R1
        lossmax = max((K - a - b for (a, b), c in R1.items() if c != 0), default=None)
        print('   %-45s R(0) = 0: %s;  dR/de in the certified kernel: %s;  antisymmetric: %s; symmetric: %s; deepest nonzero loss: %s'
              % (name, zero0, inK, anti, symm, lossmax))
        if name.startswith('A  '):
            if len(ker) == 2 and mode == 'real' and not inK:
                status = max(status, 2)
            if mode == 'tamper' and inK:
                status = 1
            if K == 3 or (len(ker) == 2):
                lead = R1.get((K, 0), 0)
                if lead != 0:
                    show('A / (X^K coefficient)', {m: c / lead for m, c in R1.items()}, K)
    if len(ker) == 2:
        lead = S[(K, 0)]
        show('S = Sol-3 table charge, monic', {m: c / lead for m, c in S.items()}, K)
        # antisymmetric element of the kernel, if any: solve (x1 e1 + x2 e2) + swap(...) = 0
        M = sp.Matrix([[e_.get(m, 0) + swap(e_).get(m, 0) for e_ in ker] for m in mons])
        ns = M.nullspace()
        print('   antisymmetric elements of the kernel: dimension %d' % len(ns))
        if ns:
            v = ns[0]
            E = {m: sum(v[i] * ker[i].get(m, 0) for i in range(2)) for m in mons}
            lead = E.get((K, 0), 0)
            if lead != 0:
                show('certified antisymmetric class / (X^K coefficient)', {m: c / lead for m, c in E.items() if c != 0}, K)
print('\nRESULT (%s):' % mode, {0: 'H-A holds at the two-dimensional kernels; controls as expected', 2: 'H-A FAILS as sealed', 1: 'control failure'}[status])
sys.exit(status)
