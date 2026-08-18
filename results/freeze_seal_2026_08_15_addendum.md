# Addendum to freeze_seal_2026_08_15.json (2026-08-18)

Two provenance fields of the seal are imprecise; per the
records-tell-the-truth rule the seal itself is NOT edited (it is a
stamped artifact); this addendum corrects the reading:

1. `launch_commit` (5988184) is the commit of the seal's RESTAMP,
   not the state at which the hours-class checks were launched.
   The checks ran at the freeze state recorded in the seal's own
   evidence entries; read `launch_commit` as
   `certificate_generation_commit` (the merged-record precedent),
   with the underlying test runs launched at the pre-restamp
   state (efc34f8 and its parent tree).
2. `worktree_note` says the tracked dirty file was "this seal
   file itself"; the recorded dirty path is actually
   results/vacancy_final_modp_smoke.json (the smoke evidence
   committed alongside).  The explanation of WHY the tree was
   dirty-by-construction stands; the file it names was wrong.

Found by the Codex release audit, 2026-08-18.
