# C₄ densest lattice packing — computational campaign

Goal: prove that the densest lattice packing of the 4-dimensional cross-polytope has determinant
`D* = 0.80822000090353537…` (density `0.82485791730…`), attained only by `Λ_ζ` and its signed-permutation images.

* Start here: `CLAUDE.md` (session brief) → `docs/STATUS.md` → `docs/PLAN.md` → `docs/DEAD_ENDS.md`.
* Regression: `python3 verify/run_all.py` (prints `ALL PASS`).
* Superseded material (do not use for proofs): `PROOF-COMPLETE.md`, `crosspolytope.tex`, `WRITEUP-SKELETON.md`,
  `dim4-*.tex`, `bnb.py`, `typebnb*.py`, `critdet.py`, `contactbf.py`, `descent_cluster.py`, `gap2_*`
  (all assert or chase the false "L₄ is optimal" / "exclude types 2–6" statements; see DEAD_ENDS A1–A4, B1–B3).
