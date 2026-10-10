# Can the Kohn–Sham scheme be modified to restore representability for OpenAI's three-electron counterexample?

Research report, 9 October 2026. Primary source: OpenAI, *A Coulomb ground-state density without Kohn–Sham ensemble representation* (25 Sep 2026, 45 pp.), read in full from the GitHub PDF (`paper.pdf`, extracted text in `paper.txt`). Semi-analytic check: `trench_model.py`.

---

## 0. Verdict in one paragraph

**"Repair KS for this density" is a theorem, not a conjecture, and the repair is minimal: no change to the KS equations, only to the admissible potential class.** The formal potential V = Δg/(2g) built in the paper is a genuine measurable function (smooth away from the two nuclei, square-integrable near them). The paper's own ground-state representation (its Eq. 5.10) shows that −½Δ + V is a nonnegative symmetric operator on C_c^∞(R³); its Friedrichs extension H_V is self-adjoint with spectrum in {0} ∪ [d, ∞), g is its simple ground state, h is an eigenfunction at d, and the Aufbau determinant ĝ↑ĝ↓h↑ is a ground state of the three-electron operator Σ_i H_V with density n (Theorem A below). So n is pure-state, Aufbau, local-potential representable, and V is the unique such potential up to a constant within the Kato–KLMN class. What fails is only that V ∉ L^{3/2} + L^∞: along the HOMO nodal plane V carries an attractive trench of fixed integrated strength √(2d) whose width shrinks like p⁻³ and whose depth grows like p³ (a δ-sheet in the limit). The paper's theorem is therefore a statement about the dual space L^{3/2}+L^∞ of Lieb's density space, i.e. about the Lagrange multiplier, not about KS theory, energies, or the T_s + E_Hxc decomposition, all of which survive intact. The one genuinely open point is whether the form domain of H_V is the standard H¹ (equivalently whether the trench strength is uniformly bounded); the strip model says yes, the paper's estimates do not reach it (Conjecture B, missing lemmas named in §3.4).

---

## 1. What the paper establishes (verified against the text)

Everything in the task brief checks out against the PDF. Points used below, with the paper's numbering:

- **Potential class and operator convention.** V = L^{3/2}(R³;R) + L^∞(R³;R); H_s(v) is "determined by its semibounded closed quadratic form" (p. 1). This convention is what the minimal fix extends.
- **Angle construction and orbital equations.** Prop. 4.1 (existence, 0 < θ < π/2 on z > 0, local W^{2,p}, p < 6); Prop. 4.3: g = √n cos θ, h = √n sin θ ∈ H¹ ∩ H²_loc, V = Δg/(2g) = Δ√n/(2√n) − |∇θ|²/2 + d sin²θ ∈ L²_loc, (−½Δ + V)g = 0, (−½Δ + V)h = dh in distributions (Eqs. 4.7–4.8).
- **Weighted gap.** Prop. 4.6, Eq. (4.23): ½∫g²|∇φ|² ≥ d(∫g²|φ|² − |∫g²φ|²/∫g²) for all φ ∈ H¹(g²).
- **Ground-state representation.** Eq. (5.10): for f ∈ E_V := {f ∈ H¹ : ∫|V||f|² < ∞}, q_V[f] := ½∫|∇f|² + ∫V|f|² = ½∫g²|∇(f/g)|² ≥ d(‖f‖² − |⟨ĝ,f⟩|²). Proved for C_c^∞ and extended to E_V by cutoffs (pp. 30–31). Also ∫n|V| < ∞ (Eq. 5.9), q_V[ĝ] = 0, q_V[h] = d.
- **Certificate and identification.** Prop. 5.1: T[Γ] + ∫V n ≥ d for every finite-kinetic-energy ensemble with density n, equality for the determinant. Prop. 5.2: any v ∈ V representing n equals V + const; the only property of v used is the form bound (5.13), ∫|v||f|² ≤ δ‖∇f‖²/2 + C_δ‖f‖², plus v ∈ L¹_loc.
- **Obstruction.** Lemma 6.1 (translated weak-* limits of the bounded part), Prop. 6.2 (every such limit is the constant c = d + k²/2), Lemma 6.3 (uniform exterior form bound), Thm. 6.4 (Agmon weight ⇒ e^{a|x|}g ∈ H¹ for some a > k; the H¹ trace on the plane contradicts √n(p,0) ≥ cA(p)p⁻⁴). Thm. 6.5: the whole family. Section 7 fixes D = 10²⁰⁰ and Z by the non-numerical formula (7.10).
- **Functional consequence.** Cor. 8.1: ∂_V T_ens[n] = ∅ and no w ∈ V represents the one-sided variation of E_Hxc; −v_ext ∈ ∂_V F_ens[n], so the interacting search is attained at n.
- **Explicit disclaimers** (p. 2): no assertion about distributional, spin-dependent, or nonlocal potentials.

**Provenance caveat.** The preprint belongs to the batch of several hundred AI-generated mathematics preprints OpenAI posted to GitHub (Nature news, 6 Oct 2026; coverage at news.lavx.hu). It is unrefereed and not on arXiv. I verified the internal logic of Sections 4–6 and 8 at the level of reading; I did **not** independently verify the Coulomb channel estimates of Section 3 (Lemmas 3.2–3.4, Prop. 3.5), which supply the on-plane tail √n(p,0) ≍ A(p)p⁻⁴ and the translation limits (3.6)–(3.7). Everything below is conditional on those estimates, exactly as the paper's Theorem 6.4 is. No published response, comment, or follow-up to this preprint was found as of 9 Oct 2026.

---

## 2. Novelty assessment

**Has anyone posed "repair KS theory for a concrete Coulomb non-representable density"?** No. Before 25 Sep 2026 there was no concrete Coulomb ground-state density proven non-representable in L^{3/2}+L^∞ (the paper's §1.1 confirms: Wagner et al. 2014 is a prescribed non-Coulomb density; Trushin–Erhard–Görling 2024 is numerical). The problem is two weeks old and nothing addressing it exists in print.

**Does an existing scheme already resolve it?** The general problem "enlarge the potential class or change the reference system so that every density is representable" is well studied. Each known scheme, checked against this case:

| Scheme | Source (checked) | Covers this case? | Why / why not |
|---|---|---|---|
| Distributional potentials, 1D torus | Sutter, Penz, Ruggenthaler, van Leeuwen, Giesbertz, J. Phys. A 57, 475202 (2024), arXiv:2312.07225 | **No** (does not apply) | Every positive H¹ density on the 1D torus is ensemble v-representable with a potential in H⁻¹ giving a self-adjoint Hamiltonian. Needs H¹ ⊂ L^∞ (1D) and a compact domain. H⁻¹(R³) potentials are not form-bounded in 3D, so the 1D dual is the wrong class in 3D. The paper disclaims it explicitly. |
| Distributional potentials, 1D interval; stability | Corso, Lett. Math. Phys. 116, 27 (2026), arXiv:2504.05501; Corso–Laestadius arXiv:2512.04726 | **No** (1D only) | Full KS theory for spinless fermions in 1D with distributional potentials; Lipschitz density-to-potential map. Same dimensional restriction. |
| Finite-temperature 1D torus | Sutter et al., J. Phys. A (2026), "v-representability on a one-dimensional torus at elevated temperatures" | **No** (1D) | Gâteaux differentiability of the thermal functional on H¹ densities. 1D only. |
| Moreau–Yosida regularisation | Kvaal, Ekström, Teale, Helgaker, JCP 140, 18A518 (2014), arXiv:1312.3734 (finite box) | **Regularises, does not represent** | For every ε > 0 the proximal density n_ε is representable by v_ε in the chosen dual; n_ε → n; since ∂F(n) = ∅ in that dual, ‖v_ε‖ → ∞. See §5(d). |
| Banach/reflexive-space KS iteration | Laestadius, Penz, Tellgren, Ruggenthaler, Kvaal, Helgaker, JCP 149, 164103 (2018), arXiv:1804.08793; Penz et al. PRL 123, 037401 (2019) | **Same as MY** | Regularised iteration converges weakly to the MY solution, not to a representing potential of n. |
| Coarse-grained representability | Lammert, JCP 125, 074114 (2006), arXiv:cond-mat/0702291; PRA 82, 012109 (2010) | **Resolves the coarse-grained problem only** | Every strictly positive cell-averaged density is uniquely ensemble representable by a cell-wise constant potential. Nothing about pointwise representability; Lammert notes the fine-grained limit need not exist, which this example instantiates. |
| Approximate / excited-state representability | Garrigue, Commun. Math. Phys. 395, 1025 (2022), arXiv:2101.01127, Cor. 3.6 | **Approximate only** | For any ρ with √ρ ∈ H¹₀(Ω) and any k there are compactly supported bounded v_n whose k-th bound *mixed* states have densities → ρ in L¹∩L^q. Exact representability is not claimed in 3D; Garrigue conjectures pure-state representable densities are not dense in 3D. |
| Non-Aufbau inversions | Trushin, Erhard, Görling, PRA 110, L020802 (2024); Erhard, Trushin, Görling, JCP 162, 034116 (2025) | **Numerical, different mechanism** | Open-shell p-block atoms invert only with Aufbau-violating occupations. The present example is closed-shell-like (b↑b↓e↑) and its mechanism is the nodal plane, not level ordering. §5(b): for even potentials the obstruction extends to non-Aufbau determinants. |
| Nodal-plane asymptotics | Gori-Giorgi, Gál, Baerends, Mol. Phys. 114, 1086 (2016), arXiv:1109.2204; Gori-Giorgi, Baerends, EPJB 91, 160 (2018), arXiv:1804.00475 | **Anticipated qualitatively** | GGB 2016 dichotomy: either the density decays with I₁ on the plane and v_s is regular, or it decays with I₀ everywhere and v_s is "strongly irregular if not divergent" at the plane. The paper proves the second branch occurs for a Coulomb system. GB 2018: the √n-equation potential diverges polynomially or exponentially on the plane. Neither paper identifies the trench structure or proposes a repaired scheme. |
| Spin-DFT; Mermin ensembles (Eschrig PRB 82, 205120 (2010)); gKS (Seidl et al. PRB 53, 3764 (1996)); 1-RDM (Gilbert 1975); Harriman construction | standard | **Not addressed to this case** | Discussed in §5. |
| Reviews | Penz et al., ACS Phys. Chem. Au 3, 334 and 492 (2023); van Leeuwen, Adv. Quantum Chem. 43, 25 (2003); Lieb, IJQC 24, 243 (1983) | — | Frame the problem; none contains a 3D Coulomb representability theorem or the Friedrichs-extension observation below. |

**Novelty conclusion.** (i) The concrete problem is new. (ii) No existing scheme resolves it. (iii) The minimal fix of §3 is new as a statement, although every ingredient is already inside the paper (its Eq. 5.10 and Prop. 4.6); the paper stops at "V ∉ V" and does not observe that V nevertheless generates a self-adjoint KS Hamiltonian.

---

## 3. The minimal fix: enlarge the potential class only

### 3.1 Statement

Notation as in the paper's §6: n, θ, g, h, V, d, k with properties (6.1)–(6.4) and the gap (4.23); these hold for the scaled density of Theorem 1.1 (and, after the dilation of Prop. 2.1, for the physical one with V_phys(r) = Z²V(Zr)).

**Theorem A (representability by the Friedrichs realisation of V).**
Let T₀ = −½Δ + V with domain C_c^∞(R³) ⊂ L²(R³).
1. T₀ is symmetric and nonnegative. Its Friedrichs extension H_V is self-adjoint, H_V ≥ 0, and its form domain Q contains E_V with closed form q̄_V[f] = ½∫g²|∇(f/g)|² there.
2. ker H_V = span{g} and σ(H_V) ∩ (0, d) = ∅; d is an eigenvalue with eigenfunction h. (No surface state of the trench lies below d or between 0 and d.)
3. H_s(V) := Σ_{i=1}^{3}(H_V ⊗ 1_{C²})_i on ∧³L²(R³;C²) is self-adjoint and bounded below, with ground energy d, and the determinant Φ = ĝ↑ ĝ↓ h↑ (ĝ = g/√2) is a ground state with n_Φ = n. Hence **n is pure-state, Aufbau, local-potential representable.**
4. If V ∈ L^{3/2}+L^∞ the operator H_V coincides with the paper's form-defined H_s(V); Theorem A extends the paper's convention rather than replacing it.
5. (Uniqueness.) Let K be the KLMN class of real u ∈ L¹_loc(R³) with ∫|u||f|² ≤ α‖∇f‖²/2 + C‖f‖² for some α < 1. If u ∈ K and H_s(u) (form sum on H¹) has a ground-state ensemble with density n, then u = V + const a.e.

**Proof.** (1) For f ∈ C_c^∞, |f|²/g is a legitimate local H¹ test for (−½Δ + V)g = 0 (g ∈ H²_loc, continuous, positive; V ∈ L²_loc), which is the paper's derivation of (5.10): ⟨f, T₀f⟩ = q_V[f] = ½∫g²|∇(f/g)|² ≥ 0. Friedrichs' theorem gives H_V. Closability is also seen directly: under f = gu the form is the weighted Dirichlet form ½∫g²|∇u|² on L²(g²dx), closable because g² is continuous and positive. E_V ⊂ Q with the same form: for f ∈ E_V the paper's cutoff argument gives χ_R f → f in H¹ and ∫|V||(1−χ_R)f|² → 0 by dominated convergence, so q_V[χ_R f − f] → 0; compactly supported H¹ functions are approximated by C_c^∞ in q_V-norm since V ∈ L²_loc and H¹ → L⁴ on compacts.
(2) g, h ∈ E_V (∫|V|(g²+h²) = ∫|V|n < ∞). The distributional equations say q_V(φ, g) = 0 and q_V(φ, h) = d⟨φ, h⟩ for φ ∈ C_c^∞; C_c^∞ is a form core, so by Kato VI.2.1(iii) g, h ∈ Dom H_V with H_V g = 0, H_V h = dh. For f ∈ C_c^∞, f/g ∈ H¹(g²), so (4.23) gives q_V[f] ≥ d(‖f‖² − |⟨ĝ,f⟩|²); both sides are continuous in form norm, hence this holds on Q, and the spectral theorem gives σ(H_V) ⊂ {0} ∪ [d,∞) with 0 simple.
(3) Second quantisation of a self-adjoint semibounded operator; the three-fermion ground energy is the sum of the three lowest min-max values counting spin, 0 + 0 + d, and Φ attains it.
(4) For V ∈ L^{3/2}+L^∞ the form q_V is closed on H¹ and C_c^∞ is a core, so the closures agree.
(5) The paper's proof of Prop. 5.2 uses v only through the form bound (5.13), which is the KLMN hypothesis, plus v ∈ L¹_loc for the division step. ∎

### 3.2 What Theorem A gives and what it does not

Kept: a local multiplicative potential (an honest function); Aufbau ground-state occupation; pure-state representation (no ensemble needed); uniqueness up to a constant within K; T_s[n] = T_ens[n] attained by Φ with T_s[n] = d − ∫Vn (paper's Eq. 5.8); hence a well-defined v_xc := V − v_ext − v_H as a function. Fixing the constant by V − c → 0 off the plane gives the standard far-field v_xc → −1/r off-plane (v_H → 3/r, v_ext → −2/r).

Not given by Theorem A alone: (a) that the form domain Q equals H¹ (equivalently V ∈ K); (b) the Banach-dual structure of Lieb's theory, since V ∉ (L¹∩L³)*. Both are addressed by Conjecture B.

### 3.3 The function space of V (proved vs conjectured)

Proved (paper + bootstrap): V ∈ L²_loc(R³), even, axisymmetric; V ∈ C^∞(R³ ∖ {±De_z}) because n is smooth and positive there (Prop. 3.1) and θ solves div(n∇θ) = −dn sin 2θ, so Δθ = −∇log n·∇θ − d sin 2θ bootstraps θ from W^{2,p}_loc to C^∞. Near the nuclei V inherits the Kato-cusp singularity of Δ√n/(2√n) (∼ −1/|x ∓ De_z|, in L²_loc). V ∉ L^{3/2}+L^∞ (Thm. 6.4).

Conjectured (strip model, §4): V is **not bounded below**; along the plane V(p,0) ≍ −(2/π)(a²/b)√(2d) p³ → −∞, and

  V(p,z) = c − √(2d) δ_w(z) + o(1),  c = d + k²/2,  w = w(p) = (b/a)p⁻³,  δ_w(z) = (2/π) w³/(z²+w²)²,  ∫δ_w = 1,

with a = |a₀'(0)| (planar slope of the direct removal amplitude, Eq. 3.26) and b = ‖T_eq‖ (planar coefficient, Eq. 3.34). So the precise class is not "bounded plus a measure on the plane" (V is a function) but

  V ∈ L^{3/2} + L^∞ + 𝒮,  𝒮 := { s : supp s ⊂ {|z| ≤ 1}, sup_{(x,y)} ∫|s(x,y,z)| dz < ∞ }.

𝒮 is infinitesimally form-bounded with respect to −Δ by the one-dimensional trace inequality sup_z|u|² ≤ η‖∂_z u‖² + η⁻¹‖u‖², so L^{3/2}+L^∞+𝒮 ⊂ K. This is the class in which "the KS equations need no modification".

**Conjecture B.** V ∈ L^{3/2} + L^∞ + 𝒮. Consequently V ∈ K, Q(H_V) = H¹, H_V is the KLMN form sum, Prop. 5.2 applies verbatim, and σ_ess(H_V) = [k²/2, ∞) (the surface-state continuum of the sheet, far above d).

### 3.4 Named unproven lemmas

- **L1 (off-plane boundedness).** sup{|V(x) − c| : |z| ≥ 1, |x| ≥ R} < ∞. Needs the vertical second derivative of log n off the plane; the paper controls only ∂_r, ∂_p and ∂_p² (Eqs. 3.3–3.4, 3.45).
- **L2 (sheet integrability).** sup_{p ≥ R} ∫_{−1}^{1} |V(p,z) − c| dz < ∞. This is the content of Conjecture B; it requires a matched inner-layer analysis of the angle equation at scale z ∼ p⁻³ with rigorous error control. Not in the paper.
- **L3 (sheet strength).** ∫_{−1}^{1}(V(p,z) − c) dz → −√(2d). Needed only for σ_ess = [k²/2,∞); it is the zero-energy bound-state condition of the limiting one-dimensional operator −½∂_z² + d − γδ(z), γ = √(2d).
- **L4 (extended subgradient = directional derivative).** For §6: that the one-sided derivative of T_ens at n along m − n equals −∫V(m − n) for all m ∈ D with ∫|V|m < ∞. The inequality "≥" is unconditional (the paper's Prop. 5.1 argument applied to arbitrary Γ with Tr γ = 3); equality is the standard convex-analysis identity, here outside a Banach setting.

Theorem A does not depend on L1–L4.

### 3.5 Where the paper's contradiction is evaded

For V of the form c − √(2d)δ_w(z): the translated weak-* limits of the bounded part are still the constant c (the sheet has vanishing width), so Prop. 6.2 is consistent; what fails is Lemma 6.3, Eq. (6.10): sup_{‖u‖_{H¹(B)} ≤ 1} |∫(V(y+x) − c)|u|²| does **not** tend to zero along the plane, because ∫δ_w|u|² → ∫_{z=0}|u|² is the trace, bounded but not small in H¹. Lemma 6.1 needs the bounded part in L^∞, which the p³-deep trench violates. The obstruction is exactly "the Lagrange multiplier has a δ-sheet component": invisible to L^∞-weak-* limits, visible to the quadratic form.

---

## 4. Semi-analytic check (strip model)

**Model.** Near the plane at large horizontal distance p, n ≈ (A(p)/p)² N(z), N(z) = a²z² + ε², ε = ε(p) = b p⁻³ (paper's Eqs. 3.5–3.6). The horizontal Laplacian acting on e^{−kp}p^m contributes −k²/2 to −½Δ, so with W(z) := V(p,z) − k²/2 the local problem is −½G'' + WG = 0, −½H'' + WH = dH, G² + H² = N, G even, H odd, with the paper's angle construction g = √N cos θ, h = √N sin θ, (Nθ')' = −dN sin 2θ, θ(0) = 0, θ(∞) = π/2, and W = (√N)''/(2√N) − θ'²/2 + d sin²θ (paper's 4.2 and 4.7). Rescaling z by a reduces to a = 1 with d → d/a², W → a²W.

**Results** (`trench_model.py`; shooting on θ'(0), DOP853, rtol 10⁻¹²; finite-difference check W = G''/(2G) to 2·10⁻⁵):

| d | ε | θ'(0)·ε | W(0)·ε | W(0)·ε / √(2d) | ∫(W − d) dz | −√(2d) | LPS (√N)''/(2√N) at 0 |
|---|---|---|---|---|---|---|---|
| 0.1 | 10⁻¹ | 1.0284 | −0.2880 | −0.644 | −0.4651 | −0.4472 | +50 |
| 0.1 | 10⁻² | 1.0028 | −0.2850 | −0.637 | −0.4491 | −0.4472 | +5·10³ |
| 0.1 | 10⁻³ | 1.0003 | −0.2847 | −0.6367 | −0.4474 | −0.4472 | +5·10⁵ |
| 0.01 | 10⁻³ | 1.0001 | −0.0900 | −0.6366 | −0.14144 | −0.14142 | +5·10⁵ |
| 0.001 | 10⁻³ | 1.00003 | −0.02847 | −0.6366 | −0.044723 | −0.044721 | +5·10⁵ |

Findings:
1. **Cancellation of the leading singularity.** (√N)''/(2√N) = +a²ε²/(2N²) is a repulsive ridge of height 1/(2ε²) ∼ p⁶ on the plane (the Gori-Giorgi–Baerends 2018 polynomial divergence of the √n-equation potential, with its sign fixed as +∞). θ'(0) = (a/ε)(1 + O(ε)) cancels it exactly at order ε⁻², leaving the KS potential at order ε⁻¹.
2. **Trench.** W − d → −√(2d)·(2/π)·ε³/(z²+ε²)² as ε → 0: depth W(0) = −(2/π)√(2d)/ε (ratio −0.6366 = −2/π to four digits), width ε, integrated strength → −√(2d) exactly, independent of a and ε. With ε = (b/a)p⁻³ this is depth ∝ p³, width ∝ p⁻³: the attractive δ-sheet. V is not bounded below.
3. **Orbital profile.** G/G(0) matches e^{−√(2d)|z|} to four digits at z = 1/√(2d) and 3/√(2d): the doubly occupied orbital is a surface state bound to the plane with binding energy d, decaying along the plane at rate k and away from it at rate √(2d). This is what lets g carry the on-plane tail A(p)p⁻⁴ while vanishing (relative to A/p) off the plane, as Prop. 6.2 requires.
4. **Semiboundedness.** −½∂_z² + W ≥ 0 by construction (G > 0 is its zero mode) although W(0) → −∞. The strength √(2d) is precisely the zero-energy bound-state condition for −½∂_z² + d − γδ(z), confirming L3 in the model.
5. **Spectrum.** In the limiting sheet model the one-particle essential spectrum starts at k²/2 ≈ 0.49 (surface continuum), far above d < 10⁻³; no state of the trench lies below d, consistent with (4.23).

No grid or basis demonstration of non-representability was attempted; by Chayes–Chayes–Ruskai every positive density on a finite lattice is representable, so none can exist.

---

## 5. Genuine modifications (needed only if one insists on L^{3/2}+L^∞)

Legend: LP = local multiplicative potential; AB = Aufbau; TS = T_s + E_Hxc decomposition with a well-defined v_xc; SC = convergent self-consistent iteration.

| Scheme | Modified problem | "Representable" means | For the OpenAI density | LP | AB | TS | SC |
|---|---|---|---|---|---|---|---|
| (0) Minimal fix (§3) | H_V = Friedrichs ext. of −½Δ+V, V ∈ K | ground ensemble of Σ_i H_V has density n | **Resolves** (Thm. A); unique in K | yes | yes | yes (v_xc a function, not in L^{3/2}+L^∞) | standard for form-bounded v; open in general |
| (a) Spin-dependent local potentials | v_↑, v_↓ ∈ V; H = Σ_i(−½Δ + v_{σ_i}) | ground ensemble with spin-summed density n | **Undecided, plausibly resolves** | yes (two) | yes | yes (v_xc^σ) | as usual |
| (b) Non-Aufbau reference, v ∈ V even | determinant of eigenfunctions of −½Δ+v, any levels | n = Σ_occ|φ_i|² | **Does not resolve** (Prop. C) | yes | no | T_s not the constrained-search value | no variational principle |
| (b') Non-Aufbau, v ∈ V not reflection-symmetric, orbitals without parity | same | same | **Undecided** | yes | no | no | no |
| (c) Finite-temperature (Mermin) ensemble | Gibbs state of H_s(v) at T > 0, v ∈ V, on R³ | thermal density = n | **Ill-posed on R³** (no trace-class Gibbs state with continuous spectrum); in a box it is a different density; continuum thermal representability open (proved only on the 1D torus) | yes | n/a (FD occupations) | yes (free-energy form) | yes (finite dim.) |
| (d) Moreau–Yosida, explicit ε | F_ε = F □ ‖·‖²/(2ε) on a Hilbert/reflexive density space X; v_ε = −∇F_ε(n) ∈ X* | n_ε := prox_ε(n) is v_ε-representable; n_ε → n | **Represents n_ε, never n**; ‖v_ε‖_{X*} → ∞ whenever X* ⊉ 𝒮 (true for X = L²: ∫V² ∼ ∫p⁴dp = ∞); v_ε → V in distributions (conjectural) | yes | yes | yes (E_Hxc,ε) | yes, weak (Laestadius 2018, Penz 2019) |
| (e) Local v ∈ V + controlled nonlocal correction near the plane | H = −½Δ + v + K, v ∈ V, K bounded, supp K ⊂ {|z| ≤ 1}² | ground ensemble density = n | **Resolves trivially or not at all**: K = (V − v)·1_trench is a multiplication operator in 𝒮 (local, form-bounded); any genuinely nonlocal K is an unsolved inverse problem; unconstrained nonlocality is trivial (Harriman/Gilbert) | no | yes | gKS-type | gKS-type |
| (f) Coarse-graining (Lammert) | cells of a partition; v constant per cell | cell averages of n_Γ equal those of n | **Resolves the coarse-grained problem** (Lammert 2006, unique potential); the fine-grained limit does not exist in V | yes (piecewise const.) | yes | yes (coarse-grained Lieb functional Gâteaux-differentiable in a box) | yes |
| (g) 1-RDM / gKS with unconstrained nonlocal potential | H = −½Δ + v + W_nl | ground-state 1-RDM or density | **Trivially resolves** (W_nl = −λ(1 − γ_Φ)); uninformative | no | yes | not a density functional | gKS |

**Prop. C (non-Aufbau with an even potential does not help).** Suppose v ∈ V is even in z and n = 2φ_a² + φ_b² with φ_a even, φ_b odd eigenfunctions of −½Δ + v at energies ε_a, ε_b (any order). Then a contradiction follows as in Thm. 6.4. Sketch: translate along the plane with the normalisation A(p_j)/p_j as in Prop. 6.2. The limits satisfy 2φ_{a,∞}² + φ_{b,∞}² = a²z²e^{−2kω·x}; φ_{a,∞} is even, vanishes on the plane, hence has zero Cauchy data there, so by unique continuation (W bounded) φ_{a,∞} ≡ 0, which replaces the Harnack step that needed positivity. Then φ_{b,∞} = ±az e^{−kω·x} forces W = ε_b + k²/2. If ε_a < ε_b, Lemma 6.3 and the Agmon argument give e^{a|x|}φ_a ∈ H¹ with a > k, contradicting the trace. If ε_a ≥ ε_b (the "hole" configuration), the natural decay rate √(2(ε_b + k²/2 − ε_a)) ≤ k of φ_a off the plane contradicts φ_a² ≤ n ∼ e^{−2kr}; this step needs a Carmona–Simon-type lower bound on eigenfunction decay (named unproven lemma L5). The paper's Prop. 5.2 handles the Aufbau case for *any* v (not necessarily even); Prop. C covers even non-Aufbau; the asymmetric non-Aufbau case (b') remains open because orbitals without parity may have translated limits c_i z e^{−kω·x} compatible with a rate-k planar trace.

**Why (a) plausibly works.** With v_↓ ≠ v_↑ the certificate of Prop. 5.1 no longer identifies the orbitals: the representing ensemble minimises T at fixed (n_↑, n_↓), not at fixed n. Route: choose a positive even g_↓ with Δg_↓/(2g_↓) ∈ L^∞ (a surface-state profile of finite width in z carrying the entire planar tail A(p)p⁻⁴ and a p⁻⁶-relative correction off-plane); set v_↓ := Δg_↓/(2g_↓). The remainder n' := n − g_↓² is positive, equals n to relative order p⁻⁶ off the plane, and decays on the plane at a rate chosen ≥ 2√(k²+2d). For n' the trace contradiction has no fuel (Eq. 6.11 needs the planar tail to beat e^{−2ap} for some a < √(k²+2d)). What remains is the inverse problem for n' with the two lowest orbitals of a single v_↑ ∈ V: an instance of the paper's own question for a non-Coulomb density with the known obstruction removed; a finite-width trench in v_↑ (allowed in L^∞, excluded in the spin-independent case only by Prop. 6.2) is available. Classification: undecided, favourable.

---

## 6. Consequences for exact and approximate functionals

**Universal planar feature of the exact v_xc.** Yes, and sharper than Gori-Giorgi–Baerends: whenever the KS HOMO has a nodal plane and the interacting density decays on it with the ground-state rate and a polynomial suppression p^{−2m} relative to the off-plane strip density (here m = 3), the exact KS potential (hence v_xc, since v_ext and v_H are regular) must carry an attractive sheet along the plane of integrated strength √(2Δε), Δε = ε_HOMO − ε_{HOMO−1}, with width ∝ p^{−m} and depth ∝ p^{m}. Its only job is to turn the doubly occupied orbital into a surface state so that it can supply an exponentially small density tail. The √n (Levy–Perdew–Sahni) potential instead carries a repulsive ridge of height ∝ p^{2m}; the two effective potentials differ in sign and order on the plane, by −|∇θ|²/2 + d sin²θ.

**Approximate functionals.** No semilocal functional should attempt it. On the plane |∇n|/n ∼ 1/w ∼ p³ and τ − τ_W ≈ ½a²(A/p)²w²/(z²+w²) is a bump of width w, so meta-GGA indicators do detect the plane, and sqrt((τ−τ_W)/n)·√(2d)/a would reproduce the depth scaling; but the prefactor contains the global quantity d (the KS gap) and the strength √(2d) is not a local functional of the density. More importantly the feature has no energetic payoff and only affects a density tail of relative size e^{−2kp}p^{−8}; a functional that reproduced it would be fitting a Lagrange-multiplier artefact. It does matter for exact-conditions work, for KS inversion algorithms (which will stall or produce wild potentials along nodal planes, a plausible contributor to the inversion failures reported by Görling's group), and for any claim that the exact v_xc of a molecule lies in L^{3/2}+L^∞.

**Energies: no effect, confirmed.** The paper gives −v_ext ∈ ∂F_ens[n], so E₀ = F[n] + ∫v_ext n with the interacting search attained; T_s[n] = T_ens[n] is attained by Φ (Prop. 5.1); E_Hxc[n] = F[n] − T_s[n] is finite and the KS energy functional T_s[m] + ∫v_ext m + E_Hxc[m] = F[m] + ∫v_ext m is minimised at m = n. Every energy quantity is exact and finite. The obstruction concerns differentiability only: no element of (L¹∩L³)* = L^∞ + L^{3/2} is a subgradient of T_ens at n (Cor. 8.1). In the extended pairing, −V is a subgradient of T_ens at n on {m ∈ D : ∫|V|m < ∞} (unconditional: the Prop. 5.1 inequality applied to arbitrary Γ with Tr γ = 3 never used n_Γ = n beyond ∫|V|n_Γ < ∞), and v_Hxc = V − v_ext represents the one-sided variation of E_Hxc along such segments modulo lemma L4. Under Conjecture B the pairing domain is all of D.

---

## 7. Research verdict

1. **Novelty.** The specific question is new (the counterexample is two weeks old, AI-generated, unrefereed, unanswered in print). The general toolbox (distributional potentials, MY, coarse-graining, non-Aufbau, spin-DFT, gKS) predates it; none of its results covers this case.
2. **The minimal fix is a theorem** (Theorem A), provable from the paper's own lemmas plus Friedrichs/KLMN machinery: n is pure-state, Aufbau, local-potential representable by the function V = Δg/(2g), unique up to a constant in the KLMN class. The KS equations are unchanged. The paper's result is accurately described as "the exact KS potential of this molecule is not in L^{3/2}+L^∞", a statement about Lieb's dual pairing.
3. **Strongest defensible statement:** "There is a real, spin-independent potential V, smooth away from the nuclei, generating a nonnegative closed form (proved) and conjecturally infinitesimally form-bounded, whose Friedrichs Hamiltonian has the Aufbau determinant with density n as a ground state; V fails to be in L^{3/2}+L^∞ only through an attractive δ-sheet of strength √(2d) along the HOMO nodal plane, with no effect on any ground-state energy." The conjectural part (V ∈ K, form domain H¹) rests on the strip model and lemmas L1–L2; the rest is proved.
4. **Genuine modifications are unnecessary** for representability. Among them, spin-DFT is the only one that plausibly restores representability *within* L^{3/2}+L^∞ (undecided); even non-Aufbau references do not help (Prop. C); finite temperature is ill-posed on R³; MY and coarse-graining regularise rather than represent; nonlocal corrections are trivial unless constrained, and then unsolved.
5. **What a follow-up paper should prove:** L1–L2 (inner-layer analysis of the angle equation along the plane), giving V ∈ L^{3/2}+L^∞+𝒮 and Q(H_V) = H¹; and L4, giving the extended Lagrange-multiplier statement. L3 and the surface-state spectral picture then follow.

---

## Sources checked

- OpenAI preprint: https://github.com/openai/math/blob/main/preprints/A-Coulomb-Ground-State-Density-without-Kohn-Sham-Ensemble-Representation-September-25-2026/paper.pdf
- Gori-Giorgi, Gál, Baerends (2016): https://arxiv.org/abs/1109.2204 · Gori-Giorgi, Baerends (2018): https://arxiv.org/abs/1804.00475
- Sutter et al. (2024): https://arxiv.org/abs/2312.07225 · finite-T torus (2026): https://jyx.jyu.fi/handle/123456789/108595
- Corso (2026): https://arxiv.org/abs/2504.05501 · Corso–Laestadius: https://arxiv.org/abs/2512.04726
- Kvaal et al. (2014): https://arxiv.org/abs/1312.3734 · Laestadius et al. (2018): https://arxiv.org/abs/1804.08793
- Lammert (2006): https://arxiv.org/abs/cond-mat/0702291
- Garrigue (2022): https://arxiv.org/abs/2101.01127
- Trushin, Erhard, Görling (2024): https://cris.fau.de/publications/327055796
- Eschrig (2010): https://arxiv.org/abs/1002.4267
- Context on the OpenAI preprint batch: https://www.nature.com/articles/d41586-026-03196-8 ; https://news.lavx.hu/article/openai-dumps-100-plus-mathematical-preprints-on-github-including-counterexamples-to-major-conjectures
