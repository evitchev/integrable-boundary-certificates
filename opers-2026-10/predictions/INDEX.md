# Registered predictions used in the paper

Each prediction file lives inside its stage directory under `tree/results/lab/large_n/` (it is not duplicated here).
"Registered" means that the SHA-256 of the file was sent to the coordinating session, and recorded, before the comparison
data were computed or opened. Times are as recorded at registration, in UTC. They are self-recorded, not third-party
timestamps. The first public commitment covering this record is the hash file of tag v1.2.1 of this repository
(2026-10-02T03:24:23Z).

| Prediction (path under tree/results/lab/large_n/) | SHA-256 | registered (UTC) | What it predicts | Status at comparison |
|---|---|---|---|---|
| vir4d_prereg/PREDICTION_VIR4d.json (seal vir4d/SEAL_VIR4d.md ab1bd1d9...) | b8446a9f28bc875c97ea80d1c6f6d5e51282b0b00b86ba5120cb2b5b9918d67c | 2026-10-01T20:04:33Z | Sol 3, spins 11 and 13, integer k | the spin-11/13 tables had not been opened by the predicting seat |
| lead_vir4d_check/VIR4R_PREDICTION_SPIN11_snapshot.json | da29b2ddac5f045072f5a64f7e45e4008a6741f60db99029d7eeeef731833fba | (snapshot of the VIR4R prediction) | Sol 3, spin 11 | registered |
| vir5/PREDICTION_VIR5a.json (seal vir5/SEAL_VIR5a.md 2efea474...) | 0c12a7e9e2108abcde2f3f651326a9eb3148d9aa82862f1ffb802b3c50df8751 | 2026-10-01T20:17:42Z | Sol 3, spins 5-13 at 9 fibres, generic t | registered by hash; the spin-11/13 tables were already open (not blind) |
| vir5c/PREDICTION_VIR5c.json | 81d4f9def62bfa064a207ff04780ebd1708ba9e9864e02fa18224abffdf4ced3 | 2026-10-01T20:56:09Z | Sol 3, spin 15 at 9 fibres | blind: the spin-15 table had not been opened by anyone |
| sol2f/PREDICTION_SOL2F.json (seal sol2f/SEAL_SOL2F.md 80274804...) | 48f1e50bfde295b6205d7c680aa7c3d8f8ab7d9d89c3a9cd5af466b7a0b0e3df | see stage NOTE | Sol 2, spin 11 | registered; the table had been opened at t = 2 only |
| vir7/PREDICTION_VIR7.json | 33af7eb7149dbcff489956c8a8fc46b6f6c9e35bde66f53b7208ba512564c229 | 2026-10-01T22:02:17Z | Sol 1, spins 5-13 at 6 fibres | spins 11/13 blind; spins 5-9 tables were open |
| exc1/PREDICTION_EXC1.json (seal exc1/SEAL_EXC1.md 7e3297d0...) | de8366b1a236ed39c8de881cb653162bf322152375e9926ca6af276165e0af3d | 2026-10-02T03:23:58Z | Sol 1, level-1 I_5 characteristic polynomials, t = 2, 9/4 | registered before the I_5 level-1 data were computed |
| exc2/PREDICTION_EXC2.json (seal exc2/SEAL_EXC2.md cfba205c...) | bc0aafcda4fe1f93b362ea9ee32e83456c9801df6b017373bdd6f6835e40bd60 | 2026-10-02T04:48:57Z | Sol 3, level-1 I_5 block, t = 2, 9/4 | blind: the I_5 block was computed afterwards |
| num1/PREDICTION_extra_M51.txt (seal num1/SEAL_NUM1.md c1418ec1...) | 0f2896c45b30279f4452df37e3c9440e9d987093d7313c150eddad57362e779a | 2026-10-01T23:58:09Z | coefficient of the extra non-WKB term at M = 51/10 | a post-hoc diagnosis, registered before the fit printed |

The checks verify these hashes where the archived scripts do (most comparison scripts assert the prediction's hash before comparing).
