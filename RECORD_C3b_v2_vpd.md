# RECORD C3b v2 — V_PD and the three-axis comparison, on Paul's ACTUAL field

**Character: correction + rerun. No new claim. Neutral ending.**

Supersedes RECORD_C3b_vpd.md, which used a constant-amplitude field missing
Paul's `psi(psi−psibar)` envelope (eq. 4.1). This rerun uses Paul's actual field
(the **original** `PaulField`, `consistency/orig/field_paul.py`, pinned to
`bf5ca06`) and the **original** C2/C3a converse-KAM/WBA code — no reimplementation.

---

## 0. LIMITATIONS (read first)

1. **Field is now Paul's spec.** `chi = iota'psi²/2 + Σ eps_mn psi(psi−psibar)
   cos(mθ−nζ)`, verbatim from the original module. B for the V_PD solver is
   `(dpsi/dzeta, dtheta/dzeta, 1)` built from `PaulField.rhs` (verified to match
   the field-line ODE to machine precision). Envelope confirmed: paul_m4
   separatrices tile `[1/8, 7/8]` exactly, Chirikov = 1.000.
2. **Geometry is still the Euclidean (ρ,θ,ζ) slab** — the same choice the
   original C2/C3a made (documented in the C2 record) and consistent with Paul's
   constant-Jacobian metric. R-weighting changes V_PD by <0.05% (RECORD_A §2.4).
3. **The V_PD solver is unchanged and previously validated** (Q1 FEM, SPD, AMG;
   max-principle + grid convergence + residuals). Only the field changed.
4. **converse-KAM / WBA are the ORIGINAL code** (`consistency/orig/`), run on
   Paul's field. So all three axes now use their own source-of-record on one
   field. No tuning.
5. **Grid data is committed this time** (`results/grid/*.npz`) to avoid the loss
   that forced the earlier reimplementation.

---

## 1. Stage 0 — why the earlier validation passed on the wrong field

The earlier Stage-1 gate passed despite the wrong field because **it probes a
single resonance, where the envelope is nearly constant**:

- Envelope `psi(psi−1)` across the core: −0.1875 (ρ=0.25) → −0.25 (ρ=0.5) →
  −0.1875 (ρ=0.75): **33% in amplitude, ~15% in island width**. Across one
  narrow island it is effectively constant, so the single-island test cannot see
  it. (Confirmed: Validation-1 transition, expressed in the *effective* local
  amplitude `eps·|psi(psi−1)|`, gives ratio 1.6–3.4 vs `eps_crit=√κ⊥/2` — the
  same as the old field; the bare-`eps` offset ×4–6 is purely the `0.25` envelope
  factor at ρ=0.5.)

- **Boundary flux (the real defect).** Angle-averaged `chi(ρ)` over the full
  domain at κ⊥=1e-6, m=12 (`results/stage0_chiprofile.npz`):

  | ρ | 0.02 | 0.05 | 0.25 | 0.50 | 0.75 | 0.95 | 0.98 |
  |---|---|---|---|---|---|---|---|
  | OLD (const amp) | 0.51 | 0.56 | 0.79 | 0.68 | 0.74 | 0.69 | 0.69 |
  | NEW (envelope)  | 0.00 | 0.00 | 0.86 | 0.85 | 0.86 | 0.28 | 0.07 |

  The OLD field has `chi ≈ 0.6` right at the Dirichlet boundaries — spurious
  parallel flux across `T=0`/`T=1`, exactly what Paul's envelope forbids. The
  NEW field's `chi → 0` at the boundaries. **Lesson: passing a validation does
  not certify the model; the validation was blind to the envelope.**

## 1b. Stage 1 — re-validation on Paul's field (the gate)

**Validation 1 (single island, m=2 n=1, resonance ρ=0.5) → PASS**

| κ⊥ | eps_trans/eps_crit (bare) | ratio in effective ε̃ | slope (√ε=0.5) |
|---|---|---|---|
| 1e-2 | 6.3 | 1.58 | — |
| 1e-4 | 6.3 | 1.58 | 0.95* |
| 1e-6 | 13.6 | 3.41 | 0.54 |

`eps_trans` scales as √κ⊥ (0.316→0.0316 over ×0.01 in κ⊥); in the effective
amplitude the transition sits within 1.6–3.4× of `eps_crit=√κ⊥/2` — a clean pass.
(*the 1e-4 slope band includes saturation.)

**Validation 2 (chaotic layer, m=12, n=2..10) → PASS**

| S | D_QL | kperp_trans (V=0.5) | ratio |
|---|---|---|---|
| 2.00 | 1.02e-3 | 1.36e-4 | 0.133 |
| 2.83 | 4.09e-3 | 4.92e-4 | 0.120 |
| 4.00 | 1.64e-2 | 1.75e-3 | 0.107 |

Transition κ⊥ ∝ quasilinear `D_QL` (constant ratio ≈ 0.12 across a 16× range) —
within one decade, same behaviour as before. **Stage-1 gate passes on Paul's
field.**

## 2. Stage 2 — four Paul fields   (filled after run)

## 3. Stage 3 — three-axis comparison, original diagnostics   (filled after run)
