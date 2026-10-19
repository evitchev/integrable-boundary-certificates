# B2 HOLD-OUT RESULT (Fable secondary seat, 2026-10-13): t = 13/5, p^2 = s = 8/9, c_{N=2} = -21/2

Run after the lead's commit of PRIMARY.md (5dec436, register item 352) and GO; PRIMARY.md sec. 9 unchanged (PRIMARY.md sha256
600bfd83...).  Started 2026-10-13T14:51:57Z (HOLDOUT_START.time), ended HOLDOUT_END.time.  Same scripts as the primary
fibres, no change; logs: towers_Z_s8-9_*, towers_K_s8-9_*, sol1_t13-5, level2_Z_t13-5, level2_K_t13-5, extras_t13-5,
products_Z_t13-5, validate2_t13-5 (full stdout).

## Registered predictions (sec. 9) versus outcome

| registered prediction | outcome | verdict |
|---|---|---|
| P-II-Z: q = 0 top pieces of the spin-3 and spin-5 projections NOT in Solution 1's commutant; pieces of weight <= 6 and constants in it | spin 3: (0,8) NO; (0,4), (0,0), (1,6), (2,4), (2,0), (4,0) yes.  Spin 5: (0,12), (0,8), (1,10), (2,8) NO; all weight <= 6 pieces and constants yes.  (Spin 4 likewise: (0,10), (1,8) NO, rest yes.) | **CONFIRMED** |
| P-II-K the same | spin 3: (0,8) NO, rest yes; spin 5: (0,12), (0,8), (1,10), (2,8) NO, rest yes | **CONFIRMED** |
| spin 2: c_2(13/5) = (t-1)/(t+3) = 2/7; (1,4) piece = -2 c_2 T_uv = -4/7 T_uv (normalisations of PRIMARY sec. 7) | c_2 = 2/7 EXACTLY; (1,4) piece = -4/7 T_uv EXACTLY | **CONFIRMED** (fourth-point test of the (1,1) form passed; c_2(t) = (t-1)/(t+3) at four fibres) |
| products reading, spin 3 (low prior): c_3 -> 32/1809, a_3 -> 2432/789 | c_3 = 16/567, a_3 = 64/21 (the decomposition T^(8) = c_3 I_3 + a_3 :T_uv T_uv: holds exactly, as it must: automatic) | **FAILED**: the (1,1)-rational forms were wrong, as the low prior anticipated; the t-dependence of (c_3, a_3) is not (1,1)-rational |

Other hold-out facts: tower dimensions as at the primary fibres (Z: 1,2,1,0,1,0,1,-,1,-,1; K: 0,2,1,0,0,0,1,-,0,-,1 at w2 = 2..12);
all charges commute with G, g', J, T; (u,v) commutant of Solution 1's vectors one-dimensional at w2 = 4, 6, 8, 10, 12 (no
extra element at this fibre); the certified I_3 converts rationally and lies in it; momentum dictionary P2 = P_X^2, Q2 = pi^2 -
rho^2 as before; the kit normalisation of I_5 is regular here and its eigenvalue equals the kit table; end-to-end eigenvalue
validation of the projection (validate2) passes for every charge of Z and K; controls (wrong dictionary s' = 65/63, perturbed
W0 momentum, perturbed Solution 1 vector) fail as required; the projected top pieces commute with S1uv, S2uv at every order and
fail against e^{-u} at the lowest orders (spin 3: order 0; spin 4: 0, 1; spin 5: 0, 1, 2); Z's I_3, I_5 differ from K's;
eigenvalue level with products: top piece = c E_3 + a E_1^2 + b E_1 + d with a != 0 (Z: c = 144/7, a = 3888/7 in the tower's
raw normalisation; canonical: c_3 = 16/567, a = 16/21).

## Statement

The PREREG's FAIL (both) branch stands at the hold-out fibre: Solution 1's tower is not the H~-vacuum projection (U(1)-coset
image as defined in the PREREG) of Zeitlin's sl^(1)(2|1) tower Z or of Kolyaskin-Litvinov's sl(2|2)^(1) tower K at spins 3, 4, 5,
at t = 2, 7/2, 11/3 and 13/5.  The agreements that do hold (weights <= 6, the spin-2 charge with c_2 = (t-1)/(t+3), the constants)
are forced: the projection commutes with S1uv, S2uv at every spin and order because P_H commutes with those two vertices, and
the commutant of that pair is one-dimensional at weights 4 and 6 (dimensions 1, 1, 2, 2, 6 at w2 = 4, 6, 8, 10, 12), so below weight 8
agreement with Solution 1 is automatic (the PREREG's automatic-agreement caveat); the only non-automatic content of Level II is
at weights >= 8, where it fails.  The failure is confined to the third screening: P_H does not commute with oint W0 because W0's
bosonised prefactor carries a J-direction term, W0 = [-(s/(4-s)) dv + (s(s-2)/(4-s)) dH~] e^{-u}, whose dH~ modes move the H~
oscillator state, so the projected charges inherit the pair's conditions and lose W0's.
