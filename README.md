# Kohn–Sham representability beyond L^{3/2}+L^∞

Working directory for the manuscript *A Kohn–Sham Hamiltonian for the OpenAI three-electron density:
representability beyond L^{3/2}+L^∞* (Cody Haycraft), a response to the OpenAI preprint
*A Coulomb ground-state density without Kohn–Sham ensemble representation* (25 Sep 2026).

## Layout
- `paper/main.tex` — manuscript body; `paper/draft.tex` is the wrapper that compiles it to `paper/draft.pdf`.
- `paper/table_strip.tex`, `paper/table_diag.tex` — Tables 1 and 2, **generated** by `python trench_model.py --table` (do not edit by hand).
- `trench_model.py` — the only numerical code (strip model; shooting with the far-field decay condition; a posteriori tail estimate, diagnostics, stability checks). Needs numpy, scipy.
- `lean/ksproof/` — Lean 4 / Mathlib project checking elementary lemmas (auxiliary; not cited in the manuscript) (`KSProof/Abstract.lean`,
  `KSProof/StripModel.lean`, `KSProof/DeltaWell.lean`; `KSProof/Axioms.lean` prints the axiom dependencies).
- `paper.pdf`, `paper.txt` — the OpenAI preprint and its extracted text (λ symbols were lost in extraction).
- `ks_repair_report.md` — the original research report (novelty survey, scheme classification).
- `tools/tectonic.exe` — portable LaTeX engine used to build the PDFs.

## Build
```
# PDFs
cd paper && ../tools/tectonic.exe draft.tex
# numerics (regenerates paper/table_strip.tex and paper/table_diag.tex)
python trench_model.py --table
# Lean (elan is installed in ~/.elan; Mathlib cache already fetched)
cd lean/ksproof && lake build && lake env lean KSProof/Axioms.lean
```

## Status (10 Oct 2026; four review rounds applied)
- Theorem A (Friedrichs realisation of V is a KS Hamiltonian; Aufbau determinant is a ground state with density n): proved from the preprint's (H1)–(H5).
- Proposition 5.4 (V bounded off the plane): proved from the preprint's channel estimates (C1)–(C3) + Gilbarg–Trudinger 9.11.
- Theorem B (form domain = H^1, KLMN, uniqueness in K): proved **conditional on Estimate 5.5** (uniform vertical L^1 bound on the strip).
- Estimate 5.5: open. Reduced (Lemma 5.6) to a horizontal-gradient bound on θ and a flux balance; Lemma 5.8 gives the sufficient flux condition in normalised form (|D_z − 4J²| ≤ K ε n²/(z²+ε²)); Remark 5.9 states the target at the plane against the exact quantity √D_z/2 (positivity of ∂²_z n on the plane for large p from the channel expansion, ε_* = ε(p)(1+O(p⁻¹)); not C(p)ε(p), which would need O(p⁻³) density asymptotics the preprint lacks) and Lemma 5.8(b) shows O(ε) relative accuracy is necessary in the model.
- Numerics (strip model only): integrated in φ = π/2 − θ (rtol 1e-12, atol 1e-30); global branch selected by the far-field decay condition at Z_bc = 12/√(2d); quantities reported on |z| ≤ 8/√(2d); the "tail" column is an a posteriori estimate, not a validated bound; the residual r of the decay condition at Z is a diagnostic (growing/decaying amplitude ratio ≈ 9r/16 at κZ = 8), not an error bound; r ≤ 4e-5 in all runs. ∫|W−d| is bounded in the computed range and approaches √(2d) from above (0.4658 at ε = 0.1, d = 0.1); the inner-profile and surface-state deviations are O(ε) (Table 2).
- Conjectures 5.10 (strip bound), 5.11 (sheet asymptotics), 5.12 (essential spectrum) stated separately; the measure limit (W−d)dz ⇀ −√(2d)δ₀ and the rescaled inner profile are distinguished and both unproved.
- Lean: `lean/ksproof` checks 17 elementary lemmas (density extension, AM–GM, Riccati, strip-model derivative system and flux identity, delta-well jump). Kept as an auxiliary check only; the manuscript no longer mentions it (removed 10 Oct 2026, since it verifies none of the analytic content).
