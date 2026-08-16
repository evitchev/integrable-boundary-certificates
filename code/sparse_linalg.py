"""Sparse fraction-free exact linear algebra for the big sectors.

The dense Fraction RREF in Sector/InvSector/MixedSector suffers
catastrophic rational-coefficient growth at dimension ~10^3 (time on
the invariant side, tens of GB on the mixed side).  Here rows are
sparse integer dicts, every operation is an integer cross-multiply
followed by content (gcd) stripping, and pivots are chosen at each
vector's maximal index -- entries stay small and zeros stay unstored.

Provides drop-in kernel computations equivalent to
genuine_kernel_inv / genuine_kernel_mixed (validated by
sparse_validate.py): vacancy dims and genuine-charge spans agree with
the dense engines up to per-vector scaling.
"""

import resource
import time
from fractions import Fraction as F
from math import gcd


class Heartbeat:
    """Rate-limited progress printer for multi-day kernel runs.
    Opt-in: constructed with label=None it is a no-op, so default
    behavior of every caller is unchanged.  Prints are single flushed
    lines tagged with the label, elapsed seconds, and peak RSS, so
    pool workers interleave legibly in a shared log."""

    def __init__(self, label, interval_s=600):
        self.label = label
        self.interval = interval_s
        self.t0 = self.last = time.time()

    def __call__(self, msg, force=False):
        if self.label is None:
            return
        now = time.time()
        if force or now - self.last >= self.interval:
            rss_mb = (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                      // 1024)
            print(f"[hb {self.label} +{int(now - self.t0)}s "
                  f"rss={rss_mb}M] {msg}", flush=True)
            self.last = now


def _strip(v):
    g = 0
    for c in v.values():
        g = gcd(g, c)
        if g == 1:
            return v
    if g > 1:
        return {k: c // g for k, c in v.items()}
    return v


class SparsePivots:
    """Integer sparse row-echelon structure keyed by max-index pivots."""

    def __init__(self):
        self.rows = {}          # pivot index -> sparse integer row dict

    def reduce(self, v):
        """Gcd-stripped remainder of v modulo the pivot span."""
        return self.reduce_scaled(v)[0]

    def reduce_scaled(self, v):
        """Returns (r, c) with r = c * (exact reduction of v modulo the
        pivot span): every cross-multiply scales by its pivot
        coefficient and every gcd strip divides, and c tracks the net
        factor so kernel coefficients can be unscaled exactly."""
        v = dict(v)
        c = F(1)
        while v:
            lead = max(v)
            if lead not in self.rows:
                g = 0
                for x in v.values():
                    g = gcd(g, x)
                    if g == 1:
                        break
                if g > 1:
                    v = {k: x // g for k, x in v.items()}
                    c = c / g
                return v, c
            row = self.rows[lead]
            a, b = row[lead], v[lead]
            new = {}
            for k, x in v.items():
                new[k] = a * x
            for k, x in row.items():
                new[k] = new.get(k, 0) - b * x
            v = {k: x for k, x in new.items() if x}
            c = c * a
            g = 0
            for x in v.values():
                g = gcd(g, x)
                if g == 1:
                    break
            if g > 1:
                v = {k: x // g for k, x in v.items()}
                c = c / g
        return v, c

    def insert(self, v):
        """Insert v; returns True if it added a new pivot."""
        r = self.reduce(v)
        if not r:
            return False
        self.rows[max(r)] = r
        return True


try:
    from math import lcm
except ImportError:  # math.lcm needs Python >= 3.9
    from math import gcd as _gcd_

    def lcm(a, b):
        return a // _gcd_(a, b) * b if a and b else 0


def _to_int_vec(cur, index):
    """Fraction current dict -> (sparse int dict, scale L) with
    intvec = L * vec."""
    L = 1
    for c in cur.values():
        L = lcm(L, F(c).denominator)
    return {index[m]: int(F(c) * L) for m, c in cur.items() if c}, L


class ModPivots:
    """GF(p) sparse row echelon, max-index pivots, rows normalized to
    leading coefficient 1.  Machine-word arithmetic: no coefficient
    growth, which removes the fraction-free echelon's tail blowup
    (the wall diagnosed on the sigma>=18 vacancy sectors)."""

    def __init__(self, p):
        self.p = p
        self.rows = {}

    def reduce(self, v):
        p = self.p
        while v:
            lead = max(v)
            row = self.rows.get(lead)
            if row is None:
                return v
            f = v[lead]
            new = dict(v)
            for k, c in row.items():
                x = (new.get(k, 0) - f * c) % p
                if x:
                    new[k] = x
                else:
                    new.pop(k, None)
            v = new
        return v

    def insert(self, v):
        r = self.reduce(v)
        if not r:
            return False
        lead = max(r)
        inv = pow(r[lead], self.p - 2, self.p)
        self.rows[lead] = {k: c * inv % self.p for k, c in r.items()}
        return True


def _modp_vec(cur, index, p):
    """Fraction current dict -> sparse GF(p) dict.  Raises if a
    denominator vanishes mod p (choose another prime)."""
    out = {}
    for m, c in cur.items():
        c = F(c)
        den = c.denominator % p
        if den == 0:
            raise ValueError("denominator divisible by p")
        v = c.numerator % p * pow(den, p - 2, p) % p
        if v:
            out[index[m]] = v
    return out


MODP_PRIME = (1 << 61) - 1  # Mersenne prime


def genuine_kernel_modp(P_low, cand_basis, Nval, residue_cur_fn,
                        gen_basis_fn, d_mono_fn, spin_fn,
                        prime=MODP_PRIME, progress=None,
                        checkpoint_path=None, checkpoint_every=1500):
    """Mod-p vacancy CERTIFIER for the genuine-kernel computation.

    Returns (tags, drank, vacant_certified) where tags = dim of the
    raw kernel mod p (candidate combinations whose residue lies in the
    mod-p d-span of the target sector) and drank = rank mod p of the
    candidate-sector d-matrix (the d-exact combinations).

    Epistemics (the directions matter):
      - solution spaces only GROW mod p:      tags >= dim W_Q;
      - matrix rank only DROPS mod p:         drank <= dim D_Q;
      - d-exact combinations are always raw:  D_Q <= W_Q.
    Hence tags == drank forces dim W_Q = dim D_Q, i.e. genuine kernel
    ZERO over Q -- a rigorous vacancy certificate from a single prime.
    Otherwise the verdict is INCONCLUSIVE with the upper bound
    genuine_Q <= tags - drank (equal to the exact dimension for all
    but finitely many primes, but not certified by one).

    Residues are computed in exact rational arithmetic (that part was
    never the bottleneck) and reduced mod p per vector; only the
    echelon runs over GF(p).

    Checkpointing (added after the sigma=20 trio was killed twice --
    at 90.5% by a session teardown, at 64.2% by a machine reboot):
    with checkpoint_path set, the candidate-sweep state (pivot rows,
    tags, next index) is pickled atomically every checkpoint_every
    candidates and on sweep completion; a matching checkpoint
    (validated against prime and candidate count) resumes the sweep
    in place, and the file is removed after the verdict returns.
    A reboot now costs at most one checkpoint interval."""
    hb = Heartbeat(progress)
    s_low = spin_fn(next(iter(P_low)))
    s_cand = spin_fn(cand_basis[0])
    nc = len(cand_basis)
    hb(f"start (mod p): {nc} candidates at spin {s_cand}", force=True)
    tgt_basis = gen_basis_fn(s_low + s_cand - 1)
    tgt_index = {m: i for i, m in enumerate(tgt_basis)}
    lower = gen_basis_fn(s_low + s_cand - 2)
    hb(f"target dim {len(tgt_basis)}; {len(lower)} d-columns",
       force=True)
    apiv = ModPivots(prime)
    for i, b in enumerate(lower):
        hb(f"d-columns {i}/{len(lower)}, rows {len(apiv.rows)}")
        col = {}
        for mono, c in d_mono_fn(b).items():
            if mono in tgt_index:
                k = tgt_index[mono] + nc
                col[k] = (col.get(k, 0) + c) % prime
        col = {k: c for k, c in col.items() if c}
        if col:
            apiv.insert(col)
    hb(f"d-span rows {len(apiv.rows)}; candidate sweep", force=True)
    tags = 0
    start_j = 0
    if checkpoint_path:
        import os
        import pickle
        if os.path.exists(checkpoint_path):
            with open(checkpoint_path, "rb") as fh:
                ck = pickle.load(fh)
            if ck["prime"] == prime and ck["nc"] == nc:
                apiv.rows = ck["rows"]
                tags = ck["tags"]
                start_j = ck["next_j"]
                hb(f"RESUMED from checkpoint at candidate {start_j}/"
                   f"{nc}, rows {len(apiv.rows)}, tags {tags}",
                   force=True)
            else:
                hb(f"checkpoint mismatch (prime/nc) -- ignored",
                   force=True)

        def _save(next_j):
            tmp = checkpoint_path + ".tmp"
            with open(tmp, "wb") as fh:
                pickle.dump({"prime": prime, "nc": nc, "tags": tags,
                             "next_j": next_j, "rows": apiv.rows},
                            fh, protocol=4)
            os.replace(tmp, checkpoint_path)
            hb(f"checkpoint written at candidate {next_j}/{nc}",
               force=True)

    for j, b in enumerate(cand_basis):
        if j < start_j:
            continue
        if (checkpoint_path and j > start_j
                and (j - start_j) % checkpoint_every == 0):
            _save(j)
        hb(f"candidate {j}/{nc}, rows {len(apiv.rows)}, tags {tags}")
        cur = residue_cur_fn(P_low, {b: F(1)}, Nval)
        aug = {k + nc: c for k, c in
               _modp_vec(cur, tgt_index, prime).items()}
        aug[j] = 1
        r = apiv.reduce(aug)
        if r and max(r) < nc:
            tags += 1
        elif r:
            lead = max(r)
            inv = pow(r[lead], prime - 2, prime)
            apiv.rows[lead] = {k: c * inv % prime for k, c in r.items()}
    if checkpoint_path:
        _save(nc)
    # rank mod p of the candidate-sector d-matrix
    cs_index = {m: i for i, m in enumerate(gen_basis_fn(s_cand))}
    span = ModPivots(prime)
    drank = 0
    for b in gen_basis_fn(s_cand - 1):
        col = {}
        for mono, c in d_mono_fn(b).items():
            col[cs_index[mono]] = (col.get(cs_index[mono], 0) + c) % prime
        col = {k: c for k, c in col.items() if c}
        if col and span.insert(col):
            drank += 1
    vacant = tags == drank
    hb(f"done: tags {tags}, d-rank {drank}, "
       f"{'VACANT (certified over Q)' if vacant else 'inconclusive: genuine_Q <= ' + str(tags - drank)}",
       force=True)
    if checkpoint_path:
        import os
        if os.path.exists(checkpoint_path):
            os.remove(checkpoint_path)
    return tags, drank, vacant


def build_dpivots(lower_basis, d_mono_fn, index):
    piv = SparsePivots()
    for b in lower_basis:
        col = {}
        for mono, c in d_mono_fn(b).items():
            if mono in index:
                col[index[mono]] = col.get(index[mono], 0) + c
        col = {k: c for k, c in col.items() if c}
        if col:
            piv.insert(col)
    return piv


def genuine_kernel_sparse(P_low, cand_basis, Nval, residue_cur_fn,
                          gen_basis_fn, d_mono_fn, spin_fn,
                          progress=None):
    """Genuine charges in span(cand_basis) commuting with oint P_low.

    Single augmented echelon: the d-columns of the target sector are
    pre-inserted tagless, candidate residues follow with identity tags
    (tag indices sit BELOW the shifted value indices, so value entries
    are always eliminated first).  A candidate combination congruent to
    zero mod d then reduces to a tag-only vector -- early-stop echelon
    detects this correctly because every entry of a span member carries
    a pivot.  Genuineness likewise: candidate-sector d-columns are
    pre-inserted, and a kernel charge is genuine iff it adds a pivot
    beyond that span."""
    hb = Heartbeat(progress)
    s_low = spin_fn(next(iter(P_low)))
    s_cand = spin_fn(cand_basis[0])
    nc = len(cand_basis)
    hb(f"start: {nc} candidates at spin {s_cand}; building target "
       f"sector (spin {s_low + s_cand - 1})", force=True)
    tgt_basis = gen_basis_fn(s_low + s_cand - 1)
    tgt_index = {m: i for i, m in enumerate(tgt_basis)}
    lower = gen_basis_fn(s_low + s_cand - 2)
    hb(f"target dim {len(tgt_basis)}; pre-inserting {len(lower)} "
       f"d-columns", force=True)
    apiv = SparsePivots()
    for i, b in enumerate(lower):
        hb(f"d-columns {i}/{len(lower)}, pivot rows {len(apiv.rows)}")
        col = {}
        for mono, c in d_mono_fn(b).items():
            if mono in tgt_index:
                k = tgt_index[mono] + nc
                col[k] = col.get(k, 0) + c
        col = {k: c for k, c in col.items() if c}
        if col:
            apiv.insert(col)
    hb(f"d-span rows {len(apiv.rows)}; candidate sweep", force=True)
    kernel_tags = []
    scales = []
    for j, b in enumerate(cand_basis):
        hb(f"candidate {j}/{nc}, pivot rows {len(apiv.rows)}, "
           f"kernel tags {len(kernel_tags)}")
        cur = residue_cur_fn(P_low, {b: F(1)}, Nval)
        iv, L = _to_int_vec(cur, tgt_index)
        scales.append(L)
        aug = {k + nc: c for k, c in iv.items()}
        aug[j] = 1
        r = apiv.reduce(aug)
        if r and max(r) < nc:
            kernel_tags.append(r)
        elif r:
            apiv.rows[max(r)] = r
    hb(f"sweep done: {len(kernel_tags)} raw kernel tags; "
       f"genuineness pass", force=True)
    # genuineness: pre-load candidate-sector d-span, count new pivots
    cs_basis = gen_basis_fn(s_cand)
    cs_index = {m: i for i, m in enumerate(cs_basis)}
    span = SparsePivots()
    for b in gen_basis_fn(s_cand - 1):
        col = {}
        for mono, c in d_mono_fn(b).items():
            col[cs_index[mono]] = col.get(cs_index[mono], 0) + c
        col = {k: c for k, c in col.items() if c}
        if col:
            span.insert(col)
    genuine = []
    for t in kernel_tags:
        cur = {}
        for j, c in t.items():
            m = cand_basis[j]
            cur[m] = cur.get(m, 0) + F(c) * scales[j]
        iv, _ = _to_int_vec(cur, cs_index)
        if span.insert(iv):
            genuine.append({m: c for m, c in cur.items() if c})
    hb(f"done: {len(genuine)} genuine kernel element(s)", force=True)
    return genuine
