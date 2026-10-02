# BR2 seal: boundary interaction from the screening pairs. Opus 5 (2), 2026-10-02. Before any computation.
Literature, read before sealing (text layers of local copies: papers/arxiv-hep-th-9604044-blz2.pdf, lvz-hep-th-0312168.pdf, lz-arXiv-1208.5259.pdf; quotes by line of the pdftotext -layout output):
- BLZ II: (2.9) V_+-(u) = :e^(+-2phi(u)):. The series (2.47) "is closely related to the boundary state in so called Boundary Sine-Gordon model" (l.593-594). The boundary action (5.1) is (kappa/beta^2) Int dt cos Phi(t, 0), and footnote 8 normalises cos Phi(t) cos Phi(t') ~ 2^-1 (i(t - t'))^(-2 beta^2), i.e. boundary dimension beta^2.
- LVZ: hairpin W-currents commute with V_+- = e^(sqrt(n) X +- i sqrt(n+2) Y) ((57)-(58)). The paperclip local IM "can be defined as the system of local operators which commute with four 'screening charges'" (l.1338-1340). The boundary state obeys (I_s - Ibar_s)|B> = 0 ((9)). The paperclip itself is a BRANE (boundary constraint (2)), not a vertex perturbation.
- LZ: the pillow is a brane model (boundary constraint, l.73, 85). No vertex-operator boundary interaction located.
## Plan and sealed statements
(1) CONVENTION (Neumann doubling, fixed by BLZ (5.1) and footnote 8 against (2.9)):
    - The boundary operator e^(k.Phi_B), with Phi_B = Phi_L + Phi_R at the boundary = 2 Phi_L, corresponds to the chiral exponential e^(2k.Phi_L), with EQUAL dimension (the chiral weight under the doubled T, including the background charge).
    - Hence B_+- := e^((1/2)(+-a X_B + beta phi_B)) <-> V_+-, and B_vir := e^((p/2) phi_B) <-> e^(p phi).
(2) FIRST-ORDER CONSERVATION THEOREM (derivation, sealed):
    - For the perturbation lambda Int dt B(t) of a Neumann boundary, d/dt of I_s + Ibar_s at O(lambda) is lambda times the boundary integral of Res_(z->t) T_s(z) B(t).
    - It is a total time derivative iff the chiral residue Res T_s(z) V(w) is a total derivative. That is VIR2's screening condition.
    - So B_+-, and B_vir, conserve the cylindrical family at first order, by VIR2 (I_3, I_5 symbolic in t; I_7, I_9 at fibres).
    - Independent machine replication with my own engine (SUPER2's eng.py) at 2 rational p per solution:
      - the commutant mod d of {V_+, V_-, e^(p phi)} at weights 4, 6 and 8 (hold-out) is 1-dim;
      - the pair alone at weight 4 is 2-dim (VIR2's number);
      - so the family is singled out by the pair INSIDE Heis_X (x) Vir_N (the Virasoro screening is automatic there).
    - Prediction: replication 90%.
(3) Controls:
    - tamper a^2 -> a^2 + delta: the commutant mod d at weight 4 drops to 0. Prediction 90%.
    - N = 1 (p^2 = 2): the Sol 1/3 pair equals LVZ's hairpin (58) at n = -1, up to normalisation, with norm 1 and orthogonal; with its mirror these are LVZ's four screening charges (pencil + machine).
    - Exact point: consistency with the record's identification of the pillow exponents with the screening data (VIR2; memory 'screening charges are our exponents'). Not an independent check.
(4) Boundary flow (sealed formulas, to be machine-checked): Delta(B_+-) = a^2/2 + beta^2/2 - beta rho:
    - Sol 1: p^2/4 = (t-1)/(2(t+1));
    - Sol 2: (p^2 - 1)/2;
    - Sol 3: 2/p^2 - 1/2, dual pair 3p^2/2 - 1;
    - Delta(B_vir) = 1 (exactly marginal direction).
    - UV boundary condition: Neumann for X and phi (phi with the Vir_N background charge). The perturbation is relevant iff 0 < Delta < 1; the t-ranges are reported.
    - IR: NOT claimed.
Caps and logs: ulimit -v 16 GB; full logs.
