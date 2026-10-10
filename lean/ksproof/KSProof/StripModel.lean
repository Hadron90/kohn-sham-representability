import Mathlib

/-!
# The strip model: the flux identity (Proposition 4.1(a),(b) of the manuscript)

Data: `N z = z² + ε²`, an angle `θ` with `N θ' = J` (the vertical flux), and the flux law
`J' = -d N sin(2θ)`, which is the one-dimensional angle equation `(N θ')' = -d N sin 2θ`.
Put `G = √N cos θ` and `H = √N sin θ`.

Step 1 (`G_hasDerivAt`, `H_hasDerivAt`, `J_hasDerivAt'`, `GH_sq`): the pair `(G, H, J)`
solves the first-order system
  `G' = (z G - J H)/N`,  `H' = (z H + J G)/N`,  `J' = -2 d G H`,  `G² + H² = N`.

Step 2 (`G'_hasDerivAt`): differentiating `G' = (zG - JH)/N` once more along the system gives
  `G'' = G (ε² - J² + 2 d N H²) / N²`.

Step 3 (`flux_identity`, `W_sub_d`): since `H² = N sin² θ`, this is
  `G'' = 2 G ((ε² - J²)/(2N²) + d sin² θ)`,   i.e.   `W - d = (ε² - J²)/(2N²) - d cos² θ`
with `W := G''/(2G)`.
-/

namespace KSProof
namespace StripModel

open Real

/-- Strip-model data. -/
structure Data where
  ε : ℝ
  d : ℝ
  θ : ℝ → ℝ
  J : ℝ → ℝ
  hε : ε ≠ 0
  hθ : ∀ z, HasDerivAt θ (J z / (z ^ 2 + ε ^ 2)) z
  hJ : ∀ z, HasDerivAt J (-d * (z ^ 2 + ε ^ 2) * Real.sin (2 * θ z)) z

namespace Data

variable (S : Data)

/-- The density profile `N z = z² + ε²`. -/
def N (z : ℝ) : ℝ := z ^ 2 + S.ε ^ 2

/-- The even orbital profile `G = √N cos θ`. -/
noncomputable def G (z : ℝ) : ℝ := Real.sqrt (S.N z) * Real.cos (S.θ z)

/-- The odd orbital profile `H = √N sin θ`. -/
noncomputable def H (z : ℝ) : ℝ := Real.sqrt (S.N z) * Real.sin (S.θ z)

/-- The first derivative of `G` along the system. -/
noncomputable def G' (z : ℝ) : ℝ := (z * S.G z - S.J z * S.H z) / S.N z

/-- The second derivative of `G` along the system. -/
noncomputable def G'' (z : ℝ) : ℝ :=
  S.G z * (S.ε ^ 2 - (S.J z) ^ 2 + 2 * S.d * S.N z * (S.H z) ^ 2) / (S.N z) ^ 2

lemma N_pos (z : ℝ) : 0 < S.N z := by
  unfold N
  have hε := S.hε
  positivity

lemma N_ne (z : ℝ) : S.N z ≠ 0 := (S.N_pos z).ne'

lemma sqrtN_pos (z : ℝ) : 0 < Real.sqrt (S.N z) := Real.sqrt_pos.mpr (S.N_pos z)

lemma sqrtN_ne (z : ℝ) : Real.sqrt (S.N z) ≠ 0 := (S.sqrtN_pos z).ne'

lemma sqrtN_mul_self (z : ℝ) : Real.sqrt (S.N z) * Real.sqrt (S.N z) = S.N z :=
  Real.mul_self_sqrt (S.N_pos z).le

lemma sqrtN_sq (z : ℝ) : Real.sqrt (S.N z) ^ 2 = S.N z :=
  Real.sq_sqrt (S.N_pos z).le

lemma N_hasDerivAt (z : ℝ) : HasDerivAt S.N (2 * z) z := by
  unfold N
  have h := (hasDerivAt_pow 2 z).add_const (S.ε ^ 2)
  convert h using 1
  simp

lemma sqrtN_hasDerivAt (z : ℝ) :
    HasDerivAt (fun x => Real.sqrt (S.N x)) (z / Real.sqrt (S.N z)) z := by
  have h := (S.N_hasDerivAt z).sqrt (S.N_ne z)
  convert h using 1
  have hs := S.sqrtN_ne z
  field_simp

lemma theta_hasDerivAt (z : ℝ) : HasDerivAt S.θ (S.J z / S.N z) z := S.hθ z

/-- `G² + H² = N`. -/
theorem GH_sq (z : ℝ) : S.G z ^ 2 + S.H z ^ 2 = S.N z := by
  unfold G H
  have h := Real.sin_sq_add_cos_sq (S.θ z)
  have hsq := S.sqrtN_sq z
  calc (Real.sqrt (S.N z) * Real.cos (S.θ z)) ^ 2 + (Real.sqrt (S.N z) * Real.sin (S.θ z)) ^ 2
      = Real.sqrt (S.N z) ^ 2 * (Real.sin (S.θ z) ^ 2 + Real.cos (S.θ z) ^ 2) := by ring
    _ = S.N z := by rw [h, hsq, mul_one]

/-- The flux law in terms of `G` and `H`: `J' = -2 d G H`. -/
theorem J_hasDerivAt' (z : ℝ) : HasDerivAt S.J (-2 * S.d * S.G z * S.H z) z := by
  have h := S.hJ z
  convert h using 1
  unfold G H
  rw [Real.sin_two_mul]
  have hsq := S.sqrtN_mul_self z
  have hN : S.N z = z ^ 2 + S.ε ^ 2 := rfl
  calc -2 * S.d * (Real.sqrt (S.N z) * Real.cos (S.θ z)) * (Real.sqrt (S.N z) * Real.sin (S.θ z))
      = -S.d * (Real.sqrt (S.N z) * Real.sqrt (S.N z)) * (2 * Real.sin (S.θ z) * Real.cos (S.θ z)) := by ring
    _ = -S.d * (z ^ 2 + S.ε ^ 2) * (2 * Real.sin (S.θ z) * Real.cos (S.θ z)) := by rw [hsq, hN]

/-- First derivative of `G = √N cos θ`: `G' = (z G - J H)/N`. -/
theorem G_hasDerivAt (z : ℝ) : HasDerivAt S.G (S.G' z) z := by
  have hc : HasDerivAt (fun x => Real.cos (S.θ x)) (-Real.sin (S.θ z) * (S.J z / S.N z)) z :=
    (S.theta_hasDerivAt z).cos
  have h := (S.sqrtN_hasDerivAt z).mul hc
  unfold G
  convert h using 1
  unfold G' G H
  have hs := S.sqrtN_ne z
  have hN := S.N_ne z
  have hsq := S.sqrtN_sq z
  field_simp
  rw [hsq]
  ring

/-- First derivative of `H = √N sin θ`: `H' = (z H + J G)/N`. -/
theorem H_hasDerivAt (z : ℝ) :
    HasDerivAt S.H ((z * S.H z + S.J z * S.G z) / S.N z) z := by
  have hs' : HasDerivAt (fun x => Real.sin (S.θ x)) (Real.cos (S.θ z) * (S.J z / S.N z)) z :=
    (S.theta_hasDerivAt z).sin
  have h := (S.sqrtN_hasDerivAt z).mul hs'
  unfold H
  convert h using 1
  unfold G
  have hs := S.sqrtN_ne z
  have hN := S.N_ne z
  have hsq := S.sqrtN_sq z
  field_simp
  rw [hsq]
  ring

/-- Second derivative of `G` along the first-order system:
`G'' = G (ε² - J² + 2 d N H²)/N²`.  This is pure algebra in `G, H, J, N`. -/
theorem G'_hasDerivAt (z : ℝ) : HasDerivAt S.G' (S.G'' z) z := by
  have hnum : HasDerivAt (fun x => x * S.G x - S.J x * S.H x)
      (1 * S.G z + z * S.G' z
        - ((-2 * S.d * S.G z * S.H z) * S.H z + S.J z * ((z * S.H z + S.J z * S.G z) / S.N z))) z :=
    ((hasDerivAt_id z).mul (S.G_hasDerivAt z)).sub ((S.J_hasDerivAt' z).mul (S.H_hasDerivAt z))
  have h := hnum.div (S.N_hasDerivAt z) (S.N_ne z)
  unfold G'
  convert h using 1
  unfold G'' G'
  have hN := S.N_ne z
  have hNdef : S.N z = z ^ 2 + S.ε ^ 2 := rfl
  field_simp
  rw [hNdef]
  ring

/-- The flux identity: `G'' = 2 G ((ε² - J²)/(2N²) + d sin² θ)`. -/
theorem flux_identity (z : ℝ) :
    S.G'' z = 2 * S.G z * ((S.ε ^ 2 - (S.J z) ^ 2) / (2 * (S.N z) ^ 2) + S.d * (Real.sin (S.θ z)) ^ 2) := by
  unfold G''
  have hN := S.N_ne z
  have hH : S.H z ^ 2 = S.N z * Real.sin (S.θ z) ^ 2 := by
    unfold H
    rw [mul_pow, S.sqrtN_sq z]
  rw [hH]
  field_simp

/-- `W - d = (ε² - J²)/(2N²) - d cos² θ` for `W = G''/(2G)`, wherever `G ≠ 0`. -/
theorem W_sub_d (z : ℝ) (hG : S.G z ≠ 0) :
    S.G'' z / (2 * S.G z) - S.d
      = (S.ε ^ 2 - (S.J z) ^ 2) / (2 * (S.N z) ^ 2) - S.d * (Real.cos (S.θ z)) ^ 2 := by
  rw [S.flux_identity z]
  have hN := S.N_ne z
  have h1 : Real.sin (S.θ z) ^ 2 = 1 - Real.cos (S.θ z) ^ 2 := by
    have := Real.sin_sq_add_cos_sq (S.θ z)
    linarith
  rw [h1]
  field_simp
  ring

end Data
end StripModel
end KSProof
