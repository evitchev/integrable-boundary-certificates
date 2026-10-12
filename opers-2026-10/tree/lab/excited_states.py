"""Excited-state matrix elements of the cylindrical charges on the two-sector Fock space (2026-09-07, v2: Y oscillators).

Modes (Euclidean cylinder coordinate x, the variable of the constants C(a, b)).  X: J(x) = P + sum_{n!=0} a_n e^{-n x}, [a_n, a_m] = n delta_{n+m}.  Y: Y^i(u) has d Y^i = Q e^i + sum_{n!=0} a^i_n e^{-inu},
[a^i_n, a^j_m] = n delta^{ij} delta_{n+m}, with the zero mode Q e (e a fixed unit vector; the momentum breaks O(N) to O(N-1)).
The letters of a density monomial are d^k J (order k) and the pairs (d^a Y . d^b Y) (two Y ends with contracted indices).
The plane-normal-ordered monomial is, as in the certified vacuum recipe (code/cyl_vev.py), the sum over partial pairings
of letters with the cylinder constants C(a, b) times the cylinder-normal-ordered product of the remaining letters, whose
matrix element between oscillator states is a Wick matching:
  * an X letter pairs with a ket oscillator a_{-n} (factor n (-in)^k), with a bra oscillator a_n (factor n (in)^k), or is
    the zero mode P (k = 0 only);
  * a Y end pairs with another Y end of a letter (constant C, the pairing above), with a ket Y-oscillator end
    (n (-in)^k, index transmitted), with a bra Y-oscillator end (n (in)^k), or is the zero mode Q e^i (k = 0 only);
  * ket oscillators may pair with bra oscillators of the same mode and species (factor n; for Y ends the index is
    transmitted);
  * every oscillator end must be matched.
Index structure: pair objects are the letters' pairs, the states' Y-oscillator pairs (singlet states a^i a^i), and free
Y ends (vector states) carrying an external index or the direction e.  Following the matching through the pair objects
gives closed loops (N each), chains from zero mode to zero mode (Q^2, e.e = 1), chains from a zero mode to an external
index (Q e^i: kept as the symbol 'ext_e'), and chains between two external ends (delta).  For the sectors computed here:
  singlet states  (X oscillators; Y pairs a^i_{-n} a^i_{-m}; 'e-vector' states a^e_{-n} = e_i a^i_{-n} are singlets of
  O(N-1) and are included in the singlet block);  transverse vector states a^perp_{-n} (index orthogonal to e).
Output: the matrix on a level's singlet block and on its transverse-vector block, eigenvalues, and the commutator test.
Checks: level 0 = certified vacuum eigenvalue; [I_3, E] = 0 on every block (exact); the X-only entries agree with
excited_x.py."""
import sys, os, json
from fractions import Fraction as F
from collections import Counter
import sympy as sp
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.environ.get("IB_CODE", os.path.join(ROOT, 'code')))
sys.path.insert(0, os.environ.get("IB_LAB", os.path.join(ROOT, 'lab')))
import vev_genfun_fast as G
Pm, Qm, Nsym = sp.symbols('P Q N')

# ---- states ----
# an "end" is (species, mode, order, obj, side) ; obj groups the two ends of a Y pair (letter pair or state pair)
# state description: X oscillators tuple; Y pair oscillators list of (n, m); Y e-oscillators tuple (modes along e);
# Y transverse oscillator (single mode, external index) or None
def partitions(n, maxpart=None):
    if maxpart is None: maxpart = n
    if n == 0: return [()]
    out = []
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k): out.append((k,) + rest)
    return out

def singlet_states(level):
    """singlet (O(N-1)) states at the level: X oscillators, Y pairs (a^i_{-n} a^i_{-m}, n >= m), e-oscillators a^e_{-n}"""
    states = []
    for lx in range(level + 1):
        for xpart in partitions(lx):
            rem = level - lx
            # distribute rem over Y pair objects (each a pair (n, m) with n + m weight) and e-oscillators
            def ypairs_and_e(r):
                out = []
                for le in range(r + 1):
                    for epart in partitions(le):
                        r2 = r - le
                        # multisets of pairs with total weight r2
                        def pairsets(w, maxpair=None):
                            if w == 0: return [()]
                            res = []
                            cands = [(n, m) for n in range(1, w + 1) for m in range(1, n + 1) if n + m <= w]
                            if maxpair is not None: cands = [c for c in cands if c <= maxpair]
                            for pr in cands:
                                for rest in pairsets(w - sum(pr), pr): res.append((pr,) + rest)
                            return res
                        for ps in pairsets(r2): out.append((ps, epart))
                return out
            for ps, epart in ypairs_and_e(rem):
                states.append((xpart, ps, epart))
    return states

def transverse_states(level):
    """states with exactly one transverse Y oscillator a^perp_{-n} (external index orthogonal to e) times a singlet remainder"""
    out = []
    for n in range(1, level + 1):
        for (xpart, ps, epart) in singlet_states(level - n): out.append((xpart, ps, epart, (n,)))
    return out

def gram(state):
    xpart, ps, epart = state[:3]
    g = sp.Integer(1)
    for n, mult in Counter(xpart).items(): g *= n ** mult * sp.factorial(mult)
    for n, mult in Counter(epart).items(): g *= n ** mult * sp.factorial(mult)
    # Y pairs a^i_{-n} a^i_{-m}: <0| a^j_m a^j_n a^i_{-n} a^i_{-m} |0> = n m N + (n m if n == m) ... computed by the engine itself
    return g

# ---- the Wick engine ----
def matrix_element(mono, bra, ket, N):
    """<bra| :mono:_plane zero mode |ket> for singlet states (X part, Y pairs, e-oscillators); returns a sympy polynomial in
    P, Q (bare) and N.  Enumerates: letter pairings (C), then matchings of all ends."""
    xs, yp = mono
    # ends: list of dicts
    ends = []
    obj = 0
    for k in xs: ends.append(dict(sp='X', k=k - 1, obj=None, side='L'))
    for (k, l) in yp:
        ends.append(dict(sp='Y', k=k - 1, obj=obj, side='L')); ends.append(dict(sp='Y', k=l - 1, obj=obj, side='L')); obj += 1
    def add_state_ends(state, side):
        nonlocal obj
        xpart, ps, epart = state[:3]; perp = state[3] if len(state) > 3 else ()
        for n in xpart: ends.append(dict(sp='X', n=n, obj=None, side=side))
        for (n, m) in ps:
            ends.append(dict(sp='Y', n=n, obj=obj, side=side)); ends.append(dict(sp='Y', n=m, obj=obj, side=side)); obj += 1
        for n in epart: ends.append(dict(sp='Y', n=n, obj=('e', obj), side=side)); obj += 1
        for n in perp: ends.append(dict(sp='Y', n=n, obj=('perp', obj), side=side)); obj += 1
    add_state_ends(bra, 'B'); add_state_ends(ket, 'K')
    L = [i for i, e in enumerate(ends) if e['side'] == 'L']
    total = sp.Integer(0)
    # enumerate perfect matchings of a subset structure: every non-letter end matched; letter ends matched or zero-mode
    n_ends = len(ends)
    def weight(i, j):
        a, b_ = ends[i], ends[j]
        if a['sp'] != b_['sp']: return None
        sa, sb = a['side'], b_['side']
        if sa == 'L' and sb == 'L':
            return sp.Rational(G.Cc(a['k'], b_['k']))
        if {sa, sb} == {'L', 'K'}:
            let, osc = (a, b_) if sa == 'L' else (b_, a)
            n = osc['n']; return n * (-n) ** let["k"]
        if {sa, sb} == {'L', 'B'}:
            let, osc = (a, b_) if sa == 'L' else (b_, a)
            n = osc['n']; return n * n ** let["k"]
        if {sa, sb} == {'B', 'K'}:
            return a['n'] if a['n'] == b_['n'] else None
        return None      # K-K or B-B never pair
    def rec(unmatched, pairs):
        if not unmatched:
            yield list(pairs); return
        i = unmatched[0]; rest = unmatched[1:]
        # option: i is a letter end left as a zero mode
        if ends[i]['side'] == 'L' and ends[i]['k'] == 0:
            yield from rec(rest, pairs + [(i, None)])
        for j in rest:
            w = weight(i, j)
            if w is None or w == 0: continue
            yield from rec([r for r in rest if r != j], pairs + [(i, j, w)])
    for pairing in rec(list(range(n_ends)), []):
        val = sp.Integer(1)
        # X zero modes and Y index structure
        # build the graph on Y objects: each Y end is a node; pairings connect nodes; objects connect their two ends
        zero_x = 0; ynodes = {}
        adj = {}
        for pr in pairing:
            if pr[1] is None:
                i = pr[0]
                if ends[i]['sp'] == 'X': zero_x += 1
                else: ynodes[i] = 'zero'
            else:
                i, j, w = pr; val *= w
                if ends[i]['sp'] == 'Y': adj.setdefault(i, []).append(j); adj.setdefault(j, []).append(i); ynodes.setdefault(i, 'pair'); ynodes.setdefault(j, 'pair')
        val *= Pm ** zero_x
        # objects: connect the two ends of each Y object (letter pair or state pair); e-oscillators are objects of one end whose other "end" is the direction e (a zero-mode-like terminus)
        obj_ends = {}
        for i, e in enumerate(ends):
            if e['sp'] == 'Y' and e['obj'] is not None:
                obj_ends.setdefault(e['obj'], []).append(i)
        # traverse chains/loops over Y ends: edges = pairings (adj) and object links
        objlink = {}
        for o, es in obj_ends.items():
            if isinstance(o, tuple) and o[0] in ('e', 'perp'):
                objlink[es[0]] = o[0]
            else:
                a_, b_ = es; objlink[a_] = b_; objlink[b_] = a_
        seen = set(); ok = True
        for start in [i for i, e in enumerate(ends) if e['sp'] == 'Y']:
            if start in seen: continue
            # walk alternating: object link, pairing link
            comp = []; frontier = [start]
            while frontier:
                v = frontier.pop()
                if v in seen: continue
                seen.add(v); comp.append(v)
                lk = objlink.get(v)
                if lk not in (None, 'e', 'perp') and lk not in seen: frontier.append(lk)
                for w_ in adj.get(v, []):
                    if w_ not in seen: frontier.append(w_)
            # classify: count termini: zero-mode ends (unpaired letter Y ends) and e-termini
            termini = [v for v in comp if ynodes.get(v) == 'zero'] + [v for v in comp if objlink.get(v) in ('e', 'perp')]
            if len(termini) == 0: val *= Nsym            # closed loop
            elif len(termini) == 2:
                kinds = sorted('zero' if ynodes.get(v) == 'zero' else objlink[v] for v in termini)
                if 'perp' in kinds:
                    if kinds != ['perp', 'perp']: ok = False; break     # perp.e = 0, perp.(Q e) = 0
                else:
                    val *= Qm ** kinds.count('zero')        # zero-mode ends carry Q e; e-termini carry e; e.e = 1
            else:
                ok = False; break
        if not ok: continue
        total += val
    return sp.expand(total)

def charge_matrix(dens, N, states):
    H = sp.zeros(len(states), len(states)); Gm = sp.zeros(len(states), len(states))
    ident = {((), ()): F(1)}
    for i, bra in enumerate(states):
        for j, ket in enumerate(states):
            Gm[i, j] = matrix_element(((), ()), bra, ket, N).subs(Nsym, sp.Rational(N.numerator, N.denominator))
            tot = sp.Integer(0)
            for mono, coef in dens.items():
                tot += sp.Rational(coef.numerator, coef.denominator) * matrix_element(mono, bra, ket, N)
            H[i, j] = sp.expand(tot.subs(Nsym, sp.Rational(N.numerator, N.denominator)))
    return H, Gm

def to_article(expr):
    Pa, Qa = sp.symbols('P_art Q_art')
    e = sp.expand(expr.subs({Pm: sp.I * Pa, Qm: sp.I * Qa}, simultaneous=True))
    return sp.expand(e.subs({Pa: Pm, Qa: Qm}, simultaneous=True))

if __name__ == "__main__":
    from mixed_engine import gen_mixed_basis, cyl_P4
    N, s_ = F(-25, 2), F(-45, 2)
    rec = json.load(open(os.path.join(os.environ.get("IB_ROOT", ROOT), 'results', 'lab', 'classical', 'kernel_w6_sol3_t5_3.json')))
    BASIS = gen_mixed_basis(6)
    ker = [{BASIS[int(i)]: F(c) for i, c in w['charge_v'].items() if F(c)} for w in rec['witnesses']]
    def member(lam):
        d = {}
        for coef, v in zip((F(1), F(lam.p, lam.q)), ker):
            for m, c in v.items(): d[m] = d.get(m, F(0)) + coef * c
        return {m: c for m, c in d.items() if c}
    charges = {"I_3": cyl_P4(3, N, s_), "e5": member(sp.Rational(-1, 6)), "Y-type": member(sp.Rational(-5, 12)), "X-type": member(sp.Rational(-7, 15))}
    level = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    states = singlet_states(level); print(f"level {level}: singlet states {states}", flush=True)
    mats = {}
    for name, dens in charges.items():
        H, Gm = charge_matrix(dens, N, states)
        if level == 0:
            print(f"   {name}: level 0 = certified vacuum value: {sp.expand(to_article(H[0, 0]) - G.vev_polynomial(dens, N)) == 0}")
        M = sp.simplify(Gm.inv() * H); mats[name] = M
    for name, M in mats.items():
        comm = sp.simplify(M * mats["I_3"] - mats["I_3"] * M)
        print(f"   [{name}, I_3] = 0 on the singlet block: {comm == sp.zeros(*comm.shape)}", flush=True)
        ev = M.eigenvals()
        print(f"   {name}: singlet eigenvalues:")
        for e_, mult in ev.items(): print("      ", sp.factor(to_article(sp.expand(e_))), f"(x{mult})", flush=True)
    tstates = transverse_states(level); print(f"level {level}: transverse-vector states {tstates}", flush=True)
    tmats = {}
    for name, dens in charges.items():
        H, Gm = charge_matrix(dens, N, tstates); tmats[name] = sp.simplify(Gm.inv() * H)
    for name, M in tmats.items():
        comm = sp.simplify(M * tmats["I_3"] - tmats["I_3"] * M)
        print(f"   [{name}, I_3] = 0 on the transverse block: {comm == sp.zeros(*comm.shape)}", flush=True)
        for e_, mult in M.eigenvals().items(): print(f"   {name}: transverse eigenvalue", sp.factor(to_article(sp.expand(e_))), f"(x{mult})", flush=True)
    print("DONE")
