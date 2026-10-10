import Mathlib

/-!
# Abstract lemmas used in the operator-theoretic argument (Section 3 of the manuscript)

* `le_of_dense_of_continuous`: an inequality between continuous functions that holds on a
  dense subset holds everywhere.  This is the step by which the gap inequality (Lemma 3.6)
  and the form identity (Lemma 3.4) pass from the core `C_c^∞` to the whole form domain.
* `two_mul_le_eta`: the weighted AM–GM inequality `2ab ≤ η a² + b²/η`, the elementary
  inequality behind the one-dimensional trace bound in Lemma 5.2.
* `riccati`: for a positive twice-differentiable `G` and `u = G'/G`,
  `G''/(2G) = (u' + u²)/2`, the Riccati form of the strip potential used in Section 4.
-/

namespace KSProof

/-- Extension of a pointwise inequality from a dense subset. -/
theorem le_of_dense_of_continuous {X : Type*} [TopologicalSpace X]
    {s : Set X} (hs : Dense s) {f g : X → ℝ} (hf : Continuous f) (hg : Continuous g)
    (h : ∀ x ∈ s, g x ≤ f x) : ∀ x, g x ≤ f x := by
  have hclosed : IsClosed {x | g x ≤ f x} := isClosed_le hg hf
  have hsub : s ⊆ {x | g x ≤ f x} := fun x hx => h x hx
  have hcl : closure s ⊆ {x | g x ≤ f x} := hclosed.closure_subset_iff.mpr hsub
  intro x
  have hx : x ∈ closure s := by
    rw [hs.closure_eq]
    exact Set.mem_univ x
  exact hcl hx

/-- Weighted AM–GM: `2ab ≤ η a² + b²/η` for `η > 0`. -/
theorem two_mul_le_eta {a b η : ℝ} (hη : 0 < η) : 2 * a * b ≤ η * a ^ 2 + b ^ 2 / η := by
  have hne : η ≠ 0 := hη.ne'
  have key : η * a ^ 2 + b ^ 2 / η - 2 * a * b = (η * a - b) ^ 2 / η := by
    field_simp
    ring
  have hnn : 0 ≤ (η * a - b) ^ 2 / η := div_nonneg (sq_nonneg _) hη.le
  linarith

/-- Riccati form: if `G` is twice differentiable with `G z ≠ 0` and `u = G'/G`, then
`u' = G''/G - u²`, i.e. `G''/(2G) = (u' + u²)/2`. -/
theorem riccati {G G' : ℝ → ℝ} {G'' z : ℝ}
    (hG : ∀ x, HasDerivAt G (G' x) x) (hG' : HasDerivAt G' G'' z) (hz : G z ≠ 0) :
    HasDerivAt (fun x => G' x / G x) (G'' / G z - (G' z / G z) ^ 2) z := by
  have h := hG'.div (hG z) hz
  convert h using 1
  field_simp

end KSProof
