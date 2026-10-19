"""towers.py -- B2: commutants (mod total derivatives) of Tower Z = {S1, S2, W0} and Tower K = {S1, S2, S3, S4}
in the full two-superfield monomial space ('all' sector, both Z2 x Z2 parities), at doubled weights w2, at rational s = p^2.

Usage: towers.py <tower Z|K> <s as fraction> <w2 list comma-separated> [tag]
Writes towers_<tower>_s<..>_w<..>.json with the kernel densities (Fractions as strings) and stats; full stdout.
"""
import sys
sys.dont_write_bytecode = True
import json
import resource
import time
from fractions import Fraction as F

import gm
from gm import joint_kernel
from super2_screening import sector_candidates
import n2frame as N2


def vertices(tower, sv, perturb=None):
    S = N2.kl_vertices(sv)
    if tower == 'K':
        V = [S[1], S[2], S[3], S[4]]
    elif tower == 'Z':
        bottom, top = N2.w0_vertex(sv, 'G')
        V = [S[1], S[2], top]
    else:
        raise ValueError(tower)
    if perturb:
        # negative control: shift the momentum of the last vertex by eps along X_0
        from super2_screening import Vertex
        k, eps = perturb
        v = V[k]
        V[k] = Vertex((v.b[0] + eps, v.b[1]), v.W, label=v.label + f'+({eps},0)')
    return V


def cur_to_json(cur):
    return [[list(map(list, m)), str(c)] for m, c in sorted(cur.items())]


def cur_from_json(lst):
    return {tuple(tuple(l) for l in m): F(c) for m, c in lst}


def main():
    tower, svs, w2s = sys.argv[1], sys.argv[2], sys.argv[3]
    tag = sys.argv[4] if len(sys.argv) > 4 else ''
    perturb = None
    if tag.startswith('perturb'):
        _, k, eps = tag.split(':')
        perturb = (int(k), F(eps))
    sv = F(svs)
    eta = N2.eta_of(sv)
    V = vertices(tower, sv, perturb)
    print(f"tower {tower}, s = {sv} (c = {3 - 12 / sv}), eta = {eta}, vertices: " + ', '.join(str(v) for v in V), flush=True)
    for v in V:
        print(f"   {v.label}: b = {v.b}, W = {v.W}", flush=True)
    out = dict(tower=tower, s=str(sv), vertices=[dict(label=v.label, b=[str(x) for x in v.b], W=cur_to_json(v.W)) for v in V],
               results={})
    for w2 in [int(x) for x in w2s.split(',')]:
        t0 = time.time()
        cands = sector_candidates(w2, 'all')
        st = {}
        ker = joint_kernel(cands, V, eta, stats=st)
        dt = time.time() - t0
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
        print(f"w2 = {w2}: commutant dim {len(ker)}  (cands {st['n_cand']}, rows {st['n_rows']}, pivots {st['n_pivots']}; "
              f"{dt:.1f}s, peak RSS {rss:.0f} MB)", flush=True)
        for i, k in enumerate(ker):
            print(f"   charge {i}: {len(k)} monomials; sample: {sorted(k.items())[:4]}", flush=True)
        out['results'][w2] = dict(dim=len(ker), secs=round(dt, 2), rss_MB=round(rss), stats=st,
                                  kernel=[cur_to_json(k) for k in ker])
    fn = f"towers_{tower}_s{svs.replace('/', '-')}_w{w2s.replace(',', '_')}{('_' + tag.replace(':', '_')) if tag else ''}.json"
    json.dump(out, open(fn, 'w'), indent=1)
    print('written', fn)


if __name__ == '__main__':
    main()
