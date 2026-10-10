# KSProof — Lean 4 formalisation of the elementary lemmas

Build: `lake build` (requires the Mathlib olean cache; `lake exe cache get` if missing).
Axiom audit: `lake env lean KSProof/Axioms.lean` — every theorem depends only on `propext`, `Classical.choice`, `Quot.sound`.

| File | Statements | Manuscript |
|---|---|---|
| `KSProof/Abstract.lean` | `le_of_dense_of_continuous`, `two_mul_le_eta`, `riccati` | Lemmas 3.4, 3.6; Lemma 5.2; §4 |
| `KSProof/StripModel.lean` | `Data.G_hasDerivAt`, `H_hasDerivAt`, `J_hasDerivAt'`, `GH_sq`, `G'_hasDerivAt`, `flux_identity`, `W_sub_d` | Proposition 4.1(a),(b) |
| `KSProof/DeltaWell.lean` | `deriv_pos`, `deriv_neg`, `deriv_right_zero`, `deriv_left_zero`, `second_deriv_pos`, `zero_energy_iff`, `jump_condition` | §4.4 item 4 |

## Scope
Only the 17 elementary statements above are formalised. Mathlib (Oct 2026) has Fourier-defined Sobolev spaces
(`TemperedDistribution.MemSobolev`) and `LinearPMap` adjoints, but no weak-derivative local/weighted Sobolev
calculus, closed quadratic forms, Friedrichs extension, KLMN theorem, unbounded spectral theorem, antisymmetric
tensor products with 1-RDMs, or interior W^{2,p} elliptic estimates; Theorems A and B of the manuscript are
therefore not stated in Lean. The axiom audit covers the 17 statements only.
