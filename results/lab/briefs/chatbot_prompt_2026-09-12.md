# Prompt for web chatbots (structural/literature consultation on the cylindrical oper), 2026-09-12

I am working on the ODE/IM correspondence for a family of integrable boundary conformal field theories (free bosons on the half-line:
one boson X and an O(N) vector Y, with a boundary interaction that makes the model integrable; the "cylindrical" families of
Vitchev hep-th/0404195 (the paperclip paper is Lukyanov-Vitchev-Zamolodchikov hep-th/0312168), whose local integrals of motion I_1, I_3, I_5, ... are known exactly on a curve of
(N, s) parametrized by t). I have exact data and I am looking for the differential operator ("oper") whose spectral determinant
encodes the vacuum eigenvalues, as in Bazhanov-Lukyanov-Zamolodchikov's ODE/IM correspondence. Please reason from the facts below;
treat every "known result" you cite as a claim I will verify, and give arXiv identifiers only (no titles from memory without an
identifier; say "I don't know" rather than guess).

FACTS (all exact, verified with computer algebra; X = P^2 and Y = Q^2 are the squared X- and Y-momenta):
1. At one point of the curve, (N, s) = (28, -9), the vacuum eigenvalues through spin 11 are exactly those of a scalar second-order
   operator psi_zz = [kappa^2 u^(3/2)(1+u)^(-5/2) - X w/2 - (Y/4 + 5/16) w(1-w) + (1-w)/16] psi, w = u/(1+u), z = log u, which is the
   pillow-brane vacuum oper of Lukyanov-Zamolodchikov arXiv:1208.5259 (their eqs. (6.1)-(6.2)) with the third momentum frozen at
   P_1^2 = 1/6. At this point the charges have representatives in Vir_1(X) x Vir_28(Y) (two stress tensors, c = 1 and c = 28), the
   spin-3 density being 4[7 :T_X^2: + 6 T_X T_Y + :T_Y^2:] modulo derivatives; the first-order change of the charges along the curve
   is absorbed by a single module of the c = 1 Virasoro primary J_4 = :(dX)^4: + (5/2):(d^2 X)^2: (the "enhanced" c = 1 charge), and
   the fixed two-Virasoro embedding fails at first order.
2. At first order around that point, NO change of the scalar operator reproduces the deformation of the eigenvalues along the curve:
   not the symbol kappa^2 Lambda(1 + h sigma(w)), not any kappa^0 potential term of momentum degree <= 2 (including X^2, XY, Y^2
   terms), not any kappa^-2 Lambda^p term of momentum degree <= 2 with any power p, not an analytic kappa-odd term, nor their unions,
   for any finite combination of monomials w^i (1-w)^j with arbitrary exponents. Each exclusion is a proved integrand identity.
3. At another point, t = 3, i.e. (N, s) = (-2, -9), c = -2 (the O(N) free-fermion point), the vacuum charges are exactly two
   DECOUPLED quantum-KdV vacuum towers at c = -2 with weights Delta_X = -X/2 - 1/8 and Delta_Y = -Y/2 (so beta_X^2 = -4X,
   beta_Y^2 = 1 - 4Y in the convention Delta = (beta^2 - 1)/8); they match the radial-oscillator determinant
   Gamma(1 + beta)/Gamma((1 + beta)/2 - E/4) through spin 9; and the O(N-1)-singlet excited states decompose into c = -2 Virasoro
   modules V(Delta_X) x V(Delta_Y + n) with multiplicities sum_n m_n q^n = prod_{s >= 2} (1 - q^s)^(-floor(s/2)) (verified through
   level 5), matching Masoero-Ruzza's theorem (arXiv:2605.24563) on the completeness of ODE/IM at c = -2.
4. At the other end of the same family, t -> infinity, (N, s) = (1, 0) (the "paperclip" point), the spin-3 density is
   x_1^4 + 6 x_1^2 y_1^2 + y_1^4 - 2 x_2^2 - 2 y_2^2 = 2[(d phi_+)^4 - (d^2 phi_+)^2 + (d phi_-)^4 - (d^2 phi_-)^2] with
   phi_+- = (X +- Y)/sqrt2: again two c = -2 KdV towers, rotated by 45 degrees, and the paperclip connection coefficient factorizes
   into two radial-oscillator Gamma factors.
5. Between the two ends of that family the charges are NOT a sum of two KdV towers in any rotated/shifted momentum frame (proved
   inconsistent at spin 5), and the first-order coupling at t = 3 is neither an off-diagonal perturbation of the two oscillators nor
   a mean-field shift of their weights and central charge; all mixed X-Y coefficients vanish to first order at t = 3, the top
   mixed coefficient of the spin-(W-1) charge being C(W,2)(t-3)/(t+2W-7), so the coupling enters the leading (WKB) symbol.
6. The vacuum charges of this family are polynomials in X, Y generated (as a ring) by the spin-1 and spin-3 charges, symmetric
   under an involution exchanging shifted X and Y; the other two families need three generators (spins 1, 3, 5).

QUESTIONS (answer each separately; label each statement as "standard result", "my derivation", or "speculation"):
A. Physical identification. N = -2 bosons and c = -2: is the Y sector at N = -2 equivalent to symplectic fermions (Parisi-Sourlas),
   so that the decoupled c = -2 towers are the two copies of the free-fermion-point KdV, and is the singlet decomposition in fact 3
   the character of a known c = -2 algebra module (the singlet algebra M(1,2) / triplet W(2,3,3,3))? Give identifiers.
B. Which known ODE/IM constructions produce a spectral determinant that is a product of two c = -2 radial-oscillator determinants
   at one point of a parameter family and a scalar second-order (pillow-type) oper elsewhere, with the coupling entering the leading
   symbol at first order? Candidates I am aware of: 2 x 2 first-order linear problems of affine-Gaudin / Kondo type
   (arXiv:2003.06694, 2010.07325), higher-order scalar opers for W-algebras (arXiv:2508.20793, 2408.12917). For each candidate, state
   what its vacuum eigenvalues would look like at first order in a coupling and whether it can produce a top-degree (WKB-leading) XY
   term at first order.
C. Given fact 1 (J_4 is the tangent at the c = 1 point) and facts 3-4 (c = -2 oscillator pairs at both ends of the other family),
   what W-algebra with two generating fields (spins 2 and 4 in the momentum-degree sense) has a c = 1 and a c = -2 point of this
   kind, and what is its known oper? Does a "commutant of Vir_1 x Vir_N in the free-boson algebra" or a coset construction fit?
D. Suggest ONE concrete, checkable computation (with the expected outcome) that would discriminate between a 2 x 2 first-order
   linear problem and a scalar higher-order oper as the global object, using only vacuum eigenvalue data through spin 13 and the
   first-order data at the two anchors.

Format: numbered answers, each at most 300 words, with the label of each statement and arXiv identifiers only; no long derivations.
