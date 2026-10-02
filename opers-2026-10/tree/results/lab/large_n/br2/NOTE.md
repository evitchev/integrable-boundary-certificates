# BR2: boundary interaction from the screening pairs. Opus 5 (2), 2026-10-02. Seal 75f55ecb (sent before computing).
Exact arithmetic (own engine eng.py from SUPER2; sympy). Runs capped (ulimit -v 16 GB); full logs. Literature: local copies, quotes located in the pdftotext -layout text (line numbers there).

## (1) Convention: fixed from BLZ II
- BLZ II (5.1): A = (1/4pi beta^2) Int (Phi_t^2 - Phi_x^2) + (kappa/beta^2) Int dt cos Phi(t, 0).
- Footnote 8: cos Phi(t) cos Phi(t') ~ 2^-1 (i(t - t'))^(-2 beta^2), i.e. boundary dimension beta^2.
- (2.9): V_+-(u) = :e^(+-2 phi(u)):. (2.47) "is closely related to the boundary state in so called Boundary Sine-Gordon model" (l.593-594).
- So: NEUMANN DOUBLING. With Phi = Phi_L + Phi_R and Phi_B = 2 Phi_L at the boundary, the boundary operator e^(k.Phi_B) is the chiral exponential e^(2k.Phi_L), with EQUAL dimension (the chiral weight under the doubled T, background charge included).
- Hence B_+- = e^((1/2)(+-a X_B + beta phi_B)) <-> V_+- (item 288), and B_vir = e^((p/2) phi_B) <-> e^(p phi).
## (2) First-order conservation: the equivalence, DERIVED (lead's point (a))
- Assumptions:
  - Neumann for X and phi, phi with the Vir_N background-charge boundary term, so that T = Tbar on the boundary and the doubling continuation Tbar_s(zbar) = T_s(zbar) holds for every density built from dPhi;
  - normal ordering compatible with the doubling;
  - first order in lambda;
  - a closed boundary (radial picture, LVZ's |z| = R).
- Statement:
  - The unperturbed Neumann state obeys (I_s - Ibar_s)|N> = 0 (doubling).
  - Let |B_lambda> = P exp(lambda Oint dtau B(tau)) |N>. At O(lambda), (I_s - Ibar_s)|B> = lambda Oint dtau [I_s - Ibar_s, B(tau)] |N>.
  - By doubling, I_s - Ibar_s on boundary insertions is the contour integral of T_s around both sides of the boundary, i.e. around tau in the doubled plane. So
        (I_s - Ibar_s)|B> = 2 pi i lambda Oint dtau Res_(w->tau)[T_s(w) V(tau)] |N> + O(lambda^2).
  - This vanishes if the residue is a total tau-derivative. Conversely, for generic momenta a non-exact residue gives a nonzero operator on |N>.
  - So LVZ's condition (9) at first order <=> the chiral residue (screening) condition. On an open line the total derivative instead gives a boundary density theta_s (Ghoshal-Zamolodchikov form).
- DIRECT CHECK WITHOUT VIR2 (br2_direct.log (a)): LVZ's PRINTED paperclip density W4^(sym) (71) is conserved mod d by all FOUR hairpin exponentials e^(+-sqrt(n) X +- i sqrt(n+2) Y) ((57), (58) and the mirror). This holds at n = -18/25, -32/25 and -50/169 (rational, real-ised by X = i Xt/sqrt2), and a tamper (+1/7 on the 6n(n+2) term) fails.
  - So the residue form of (9) holds for LVZ's own data. Together with the derivation, LVZ's hairpin exponentials are first-order-conserved boundary insertions for the paperclip IM.
- OUR SOLUTIONS (independent replication of VIR2, br2_commutant.log; 2 rational fibres per solution):
  - the commutant mod d of {V_+, V_-, e^(p phi)} is 1-dim at weights 4, 6 and 8 (8 = the I_7 hold-out);
  - the pair alone gives 2 (Sols 1, 2) or 3 (Sol 3) at weight 4, matching VIR2's table;
  - TAMPER a -> a + 1/7 gives 0.
  - So B_+- conserve the family at first order, and only the triple singles it out.
## Lead's point (b): the Virasoro screening is BULK structure, not a boundary interaction
- e^(p phi) has Delta = 1 (br2_dims.log). Its charge commutes with the whole Vir_N algebra, so every density in Heis_X (x) Vir_N commutes with it automatically.
- It is not a field of the physical Vir_N CFT: it is the Coulomb-gas intertwiner encoding the phi background charge.
- So it does NOT belong to the boundary interaction. Inside the physical bulk Heis_X (x) Vir_N, the pair alone singles out the family (pair + Virasoro = 1-dim).
- Mirroring LVZ's bookkeeping:
  - Our pair (+-a, beta) is ONE hairpin's screening pair in LVZ form, with LVZ's dilaton axis = our phi (background rho) and LVZ's Y = our X.
  - For Sols 1 and 2 at EVERY t: both pairs have norm 1 and v_+.v_- = -(n+1), giving n = -p^2/2 (Sol 1) and n = -2p^2 (Sol 2). The Sol 3 pairs (norm 2/p^2, 2p^2) are not of hairpin type.
  - The mirror pair (+-a, -beta) plays the second hairpin. At N = 1 (rho = 0) it is ALSO conserved: br2_direct.log (b) shows the pair + Virasoro class at weight 4 commutes with the mirror too, so the record's I_3 there is a paperclip IM in LVZ's four-screening sense.
  - For N != 1 the mirror fails (VIR2 tamper A: rho -> -rho), and the Virasoro screening (the bulk Vir_N) takes over its role in reducing the W-commutant to a commuting family.
  - Caution: N = 1 is a kernel-jump fibre. The four-screening commutant is 2-dim at weight 4, and pair + Virasoro is 2-dim at weight 6.
## (3) Controls
- Tamper (wrong a^2): FAILS, as required.
- N = 1: the pair is LVZ's hairpin at n = -1 (norm 1, orthogonal); with its mirror these are LVZ's four screening charges, and the direct LVZ check (a) holds at rational n.
- Exact point (t = -/+1/3): Sols 1 and 2 have the SAME boundary dimension Delta = -1 there, consistent with the merger. Consistency with the pillow rests on the record's identification of the pillow exponents with the screening data. LZ describe the pillow only as a BRANE (boundary constraint, l.73, 85); no vertex-operator boundary interaction was located in LZ. So no independent pillow check.
## (4) Boundary flow (br2_dims.log)
- UV: Neumann for X and phi (phi with the Vir_N background charge), i.e. the free/Neumann boundary of Heis_X (x) Vir_N.
- Perturbation: lambda Int dt (B_+ + B_-), B_+- = e^(+-(a/2) X_B) psi_B, where psi_B = e^((beta/2) phi_B) is the boundary degenerate field of Vir_N (level 2 for Sols 1 and 3, level 3 for Sol 2 and Sol 3's dual, per VIR2).
- Dimensions:
  - Sol 1: Delta = p^2/4 = (t-1)/(2(t+1)); relevant for t > 1 or t < -3.
  - Sol 2: Delta = (p^2-1)/2 = (t-3)/(2(t+1)); relevant for t > 3 or t < -5.
  - Sol 3: Delta = 2/p^2 - 1/2 = (t+3)/(2(t-1)); relevant for t > 5 or t < -3. Its dual pair: Delta = 3p^2/2 - 1 = 2(t-2)/(t+1); relevant for 2 < t < 5.
  - Paperclip (t -> infinity): Delta -> 1/2 (Sols 1, 3). Exact point: Delta = -1 (non-unitary).
- RG direction: a relevant boundary perturbation flows away from Neumann. The IR boundary condition is NOT determined here (BR3).
- Status: first-order integrability only. Higher orders (lambda^2 counterterms, the closure of the full boundary-state equation) are NOT shown.
## Disclosures
- br2_direct.py first imported br2_commutant.py, whose unguarded main re-ran and exited before the direct checks. I fixed it with a __main__ guard and regenerated both logs.
- The derivation in (2) is analytic; the machine checks are the residue identities. I make no claim beyond first order.
