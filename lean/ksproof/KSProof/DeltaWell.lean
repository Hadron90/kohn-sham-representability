import Mathlib

/-!
# The limiting delta well

For `κ > 0` let `G z = exp (-κ |z|)`.  Then
* `deriv_pos`: on `z > 0`, `G' = -κ G`;
* `deriv_neg`: on `z < 0`, `G' = κ G`;
* `deriv_right_zero`, `deriv_left_zero`: the one-sided derivatives at `0` are `-κ` and `κ`,
  so the jump of `G'` across `0` is `-2κ`;
* `second_deriv_pos`: on `z > 0`, `G'' = κ² G`, hence `-½ G'' + d G = 0` iff `κ² = 2d`.

The jump condition for `-½ ∂² + d - γ δ(z)` reads `-½ (G'(0+) - G'(0-)) = γ G(0)`, i.e.
`γ = κ`; with `κ² = 2d` this identifies the sheet strength `γ = √(2d)` (Section 4, item 4).
-/

namespace KSProof
namespace DeltaWell

open Real Filter Topology

/-- `G z = exp (-κ |z|)`. -/
noncomputable def G (κ z : ℝ) : ℝ := Real.exp (-κ * |z|)

/-- On `z > 0` the function agrees with `exp (-κ z)` near `z`. -/
lemma G_eventuallyEq_pos (κ : ℝ) {z : ℝ} (hz : 0 < z) :
    (fun x => Real.exp (-κ * x)) =ᶠ[𝓝 z] G κ := by
  filter_upwards [lt_mem_nhds hz] with x hx
  simp [G, abs_of_pos hx]

/-- On `z < 0` the function agrees with `exp (κ z)` near `z`. -/
lemma G_eventuallyEq_neg (κ : ℝ) {z : ℝ} (hz : z < 0) :
    (fun x => Real.exp (κ * x)) =ᶠ[𝓝 z] G κ := by
  filter_upwards [gt_mem_nhds hz] with x hx
  simp [G, abs_of_neg hx]

lemma exp_neg_hasDerivAt (κ z : ℝ) :
    HasDerivAt (fun x => Real.exp (-κ * x)) (-κ * Real.exp (-κ * z)) z := by
  have h0 : HasDerivAt (fun x : ℝ => x) 1 z := hasDerivAt_id z
  have h := (h0.const_mul (-κ)).exp
  exact h.congr_deriv (by ring)

lemma exp_pos_hasDerivAt (κ z : ℝ) :
    HasDerivAt (fun x => Real.exp (κ * x)) (κ * Real.exp (κ * z)) z := by
  have h0 : HasDerivAt (fun x : ℝ => x) 1 z := hasDerivAt_id z
  have h := (h0.const_mul κ).exp
  exact h.congr_deriv (by ring)

/-- `G' = -κ G` on `z > 0`. -/
theorem deriv_pos (κ : ℝ) {z : ℝ} (hz : 0 < z) : HasDerivAt (G κ) (-κ * G κ z) z := by
  have h := exp_neg_hasDerivAt κ z
  have hG : G κ z = Real.exp (-κ * z) := by simp [G, abs_of_pos hz]
  rw [hG]
  exact h.congr_of_eventuallyEq (G_eventuallyEq_pos κ hz).symm

/-- `G' = κ G` on `z < 0`. -/
theorem deriv_neg (κ : ℝ) {z : ℝ} (hz : z < 0) : HasDerivAt (G κ) (κ * G κ z) z := by
  have h := exp_pos_hasDerivAt κ z
  have hG : G κ z = Real.exp (κ * z) := by simp [G, abs_of_neg hz]
  rw [hG]
  exact h.congr_of_eventuallyEq (G_eventuallyEq_neg κ hz).symm

/-- Right derivative at `0` is `-κ`. -/
theorem deriv_right_zero (κ : ℝ) : HasDerivWithinAt (G κ) (-κ) (Set.Ici 0) 0 := by
  have h : HasDerivWithinAt (fun x => Real.exp (-κ * x)) (-κ) (Set.Ici 0) 0 := by
    have := (exp_neg_hasDerivAt κ 0).hasDerivWithinAt (s := Set.Ici 0)
    simpa using this
  refine h.congr ?_ ?_
  · intro x hx
    simp [G, abs_of_nonneg (Set.mem_Ici.mp hx)]
  · simp [G]

/-- Left derivative at `0` is `κ`. -/
theorem deriv_left_zero (κ : ℝ) : HasDerivWithinAt (G κ) κ (Set.Iic 0) 0 := by
  have h : HasDerivWithinAt (fun x => Real.exp (κ * x)) κ (Set.Iic 0) 0 := by
    have := (exp_pos_hasDerivAt κ 0).hasDerivWithinAt (s := Set.Iic 0)
    simpa using this
  refine h.congr ?_ ?_
  · intro x hx
    simp [G, abs_of_nonpos (Set.mem_Iic.mp hx)]
  · simp [G]

/-- Second derivative on `z > 0`: the derivative of `x ↦ -κ G x` is `κ² G z`. -/
theorem second_deriv_pos (κ : ℝ) {z : ℝ} (hz : 0 < z) :
    HasDerivAt (fun x => -κ * G κ x) (κ ^ 2 * G κ z) z := by
  have h := (deriv_pos κ hz).const_mul (-κ)
  convert h using 1
  ring

/-- Zero-energy condition: `-½ G'' + d G = 0` with `G'' = κ² G` holds iff `κ² = 2d` (for `G ≠ 0`). -/
theorem zero_energy_iff (κ d g : ℝ) (hg : g ≠ 0) :
    -(1 / 2) * (κ ^ 2 * g) + d * g = 0 ↔ κ ^ 2 = 2 * d := by
  constructor
  · intro h
    have : (2 * d - κ ^ 2) * g = 0 := by linarith
    rcases mul_eq_zero.mp this with h1 | h1
    · linarith
    · exact absurd h1 hg
  · intro h
    rw [h]
    ring

/-- Jump condition for the delta well: `-½ (G'(0+) - G'(0-)) = γ G(0)` with
`G'(0+) = -κ`, `G'(0-) = κ`, `G(0) = 1` forces `γ = κ`. -/
theorem jump_condition (κ γ : ℝ) : -(1 / 2) * ((-κ) - κ) = γ * G κ 0 ↔ γ = κ := by
  simp [G]
  constructor <;> intro h <;> linarith

end DeltaWell
end KSProof
