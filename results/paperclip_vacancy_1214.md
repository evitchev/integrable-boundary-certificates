# Paperclip even-spin vacancies extended to sigma = 12, 14

Script `code/paperclip_vacancy_1214.py` (sparse stack + heartbeats),
artifact `paperclip_vacancy_1214.json`, run 2026-08-06.

| sigma | n = 7/3 | n = 13/5 | verdict |
|---|---|---|---|
| 12 | ker dim 0 | ker dim 0 | VACANT |
| 14 | ker dim 0 | ker dim 0 | VACANT |

Kernel dimensions are upper semicontinuous in n, so one generic point
proves generic vacancy; two points were run as a cross-check.  With
`vacancy.py` (sigma <= 10) this extends the paperclip even-vacancy
table to all even sigma <= 14, enlarging the certified hypothesis
range of the bootstrap lemma for the paperclip tower.

Cost note: sigma = 14 (980 candidates, target dimension 3088, 2104
d-columns) completed in ~80 s per point under `genuine_kernel_sparse`
— the sector that motivated the "needs sparse elimination" [road-map]
caveat is no longer a wall.
