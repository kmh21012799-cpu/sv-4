# RECORD C3b — V_PD implementation and the three-axis comparison

**Character: reproduction + comparison. No claim of new discovery. Neutral ending.**

Paul, Hudson & Helander (2022) computed the effective volume of parallel
diffusion `V_PD` and deferred the comparison against a converse-KAM calculation
to "future publication". C2 (converse-KAM) and C3a (WBA) established that
neither topology (converse-KAM) nor dynamics (WBA) distinguishes Paul's four
critical-overlap fields. This record adds the third axis — transport (`V_PD`) —
and puts all three on one stage: **same fields, same grid, same core.**

---

## 0. LIMITATIONS (read first)

These bound every number below. None of them is tuned away.

1. **Geometry is a reduced slab, not a torus.** Paul solves in a real toroidal
   annulus. We use the flat "screw-pinch" field-line Hamiltonian in
   `(rho, theta, zeta)` with a Euclidean metric — the geometry in which Paul's
   own boundary-layer analysis (his §5–6) is derived. Consequence: our absolute
   numbers are *not* comparable to Paul's absolute numbers; only scalings,
   orderings and qualitative structure are. The field is
   `B = (-psi_theta, psi_rho, 1)`, `psi = Psi(rho) + Σ eps_mn cos(mθ-nζ)`,
   `iota(rho) = 0.5 + rho`. Within this geometry the discretised operator is
   *exact* (no `|B|≈1` approximation): the full 3×3 tensor
   `kappa = kperp I + (kpar-kperp) BB/|B|²` is assembled directly.

2. **Nothing is tuned to Paul's values.** `iota' = 1`, `iota_0 = 0.5` are fixed
   once. Every mode amplitude follows from the Chirikov formula
   `eps = (s/4m)²/iota'`; the critical-overlap fields use `s = 1`. We did not
   adjust anything to hit `eps_crit`, the transition `kperp`, or any ordering.
   (Same discipline as B2's untuned `eps_cr`.)

3. **converse-KAM and WBA are reimplementations.** The original C2/C3a code is
   not in this repository. `vpd/diagnostics.py` contains self-contained
   field-line diagnostics: a weighted-Birkhoff `dig` (WBA) and a MacKay-style
   tangent cone-crossing converse-KAM. They reproduce the *qualitative*
   structure of the C2/C3a table (converse-KAM non-existence ≈100%, dig median
   ≈1, chaos fraction ≈90%, and the island split) but not its exact numbers
   (e.g. our `t_c` runs ~4 vs the table's ~19 — a normalisation difference).
   The comparison stands on the qualitative behaviour, not on matched digits.

4. **The Heaviside indicator is grid-sensitive (trap 2).** `V_PD` integrates a
   binary condition; `DeltaT` measures a small residual core drop. `V_PD`
   converges to ≲2%. `DeltaT` converges more slowly (~5% per refinement, still
   drifting down at our finest grid); we therefore trust `DeltaT` **orderings**
   between fields, not its absolute value. Because the per-`m` grid necessarily
   refines `Nzeta` with `m`, and finer grids give *smaller* `DeltaT`, any
   surviving "`DeltaT` larger at high `m`" result is conservative.

5. **Conditioning (trap 1).** `kpar/kperp` up to `1e6`. We assemble an SPD
   FEM stiffness matrix and solve with CG + smoothed-aggregation AMG; the AMG
   solution matches a direct sparse LU solve to 4 digits (`V_PD` 0.1381 vs
   0.1381). Residuals are `≤1e-8`. The discrete maximum principle is checked on
   every solve; interior overshoots are `<0.4%` on production grids and vanish
   under refinement (`1e-4 → 1e-8`).

---

## 1. Stage 1 — V_PD implementation and validation (the gate)

Method: bilinear/trilinear (Q1) finite elements; Fourier-spectral gradients in
the periodic angles for the `V_PD` indicator; 4th-order-equivalent convergence
checked directly. `V_PD = ∫_core Θ(kpar|∇∥T|² − kperp|∇⊥T|²) / ∫_core`.

### Validation 1 — single island chain (Paul Fig 1)   → **PASS**

`(m,n) = (2,1)`, resonance centred, `eps` scanned across ~4 decades, three
`kperp`. Prediction: `V_PD` turns positive at `eps_crit = √kperp/2`; `V_PD ∝ √eps`
at large `eps`.

| `kperp` | `eps_crit=√kperp/2` | `eps` at V_PD>0 | ratio | within 1 decade? |
|---|---|---|---|---|
| 1e-2 | 5.0e-2 | 7.0e-2 | 1.40 | yes |
| 1e-4 | 5.0e-3 | 7.0e-3 | 1.40 | yes |
| 1e-6 | 5.0e-4 | 1.4e-3 | 2.75 | yes |

- Transition sits within a factor **1.4–2.8** of `eps_crit`.
- `eps_trans` scales as **√kperp** (0.070 → 0.0070 as kperp drops 100×).
- Large-`eps` log-log slope **0.57** (predicted 0.5) at `kperp=1e-6`.
- Grid convergence: `V_PD → 0.209 ± 0.001` (<1%); overshoot `1.1e-4 → 1e-8`.

Figure: `figures/fig1_island.png`.

### Validation 2 — strongly chaotic layer (Paul Fig 3)   → **PASS**

`m=12`, resonances filling the core, Chirikov `s = 2, 2√2, 4`; `kperp` scanned
`1e-2…1e-6`. Quasilinear prediction `D_QL = Σ π (eps m)²/(m iota')`.

| `s` | `D_QL` | `kperp` at V_PD=0.5 | ratio `/D_QL` |
|---|---|---|---|
| 2.00 | 7.95e-4 | 5.05e-5 | 0.063 |
| 2.83 | 3.18e-3 | 1.94e-4 | 0.061 |
| 4.00 | 1.27e-2 | 7.42e-4 | 0.058 |

- The transition `kperp` is **proportional to `D_QL`** with a nearly constant
  ratio ≈ **0.06** across a 16× range — well inside one decade of quasilinear
  theory.
- Stronger overlap transitions at larger `kperp`; all curves rise monotonically
  from ~0 to ~1 and converge as `kperp → 0`.
- (The *numerical* map estimate `D_FL` saturates against the domain for `s=4`
  and is unreliable there; the analytic `D_QL` is the correct predictor.)

Figure: `figures/fig2_chaos.png`.

**Gate result: both analytic predictions reproduced within one order of
magnitude, with the correct scalings. Stage 1 passes; C3b proceeds.**

---

## 2. Stage 2 — the four fields (Paul Fig 5, 6)

Fields m = 4, 12, 20, 36 at critical overlap (`s=1`), core `rho ∈ [0.25,0.75]`,
`kperp = 1e-4, 1e-5, 1e-6`. Grids `Nrho=129`, `Ntheta=32` (poloidal-periodicity
reduction `theta∈[0,2π/m)`), `Nzeta` scaled to resolve the resonant `n`
(32→128). Warm-start continuation in `kperp`. Residuals `≤1e-8`, overshoot
`<1%`.

**V_PD** (parallel-diffusion volume fraction):

| `kperp` | m=4 | m=12 | m=20 | m=36 |
|---|---|---|---|---|
| 1e-4 | **0.289** | 0.103 | 0.057 | 0.038 |
| 1e-5 | **0.691** | 0.260 | 0.156 | 0.091 |
| 1e-6 | **0.884** | 0.710 | 0.465 | 0.292 |

**ΔT = ⟨T(0.75)⟩−⟨T(0.25)⟩** (temperature drop retained across the core;
large = insulates, small = leaks):

| `kperp` | m=4 | m=12 | m=20 | m=36 |
|---|---|---|---|---|
| 1e-4 | 0.455 | 0.480 | 0.485 | **0.496** |
| 1e-5 | 0.404 | 0.427 | 0.427 | **0.472** |
| 1e-6 | 0.391 | 0.399 | 0.379 | **0.440** |

Findings vs Paul §7:
- **V_PD is largest for m=4 at every `kperp`** — and monotone decreasing in `m`.
  ✔ (Paul §7)
- **ΔT is largest for m=36 at every `kperp`** — the field whose flux surfaces are
  *most* broken (highest `m`) retains the *most* temperature gradient, i.e.
  insulates best. ✔ (Paul Fig 5b, the paper's central figure). m=4 leaks most
  (smallest ΔT) at 1e-4 and 1e-5; at 1e-6 m=4 and m=20 are within the ΔT grid
  uncertainty (~5%).
- **Convergence as `kperp→0`: not seen in our range.** The V_PD spread across
  fields is 0.25 → 0.60 → 0.59 (not decreasing). Within `[1e-4,1e-6]` the four
  fields remain well separated and ordered; m=4 saturates toward 1 while high-`m`
  fields lag. The eventual convergence to `V_PD→1` lies at smaller `kperp` than
  we computed. Reported as observed, not forced.

Figure: `figures/fig3_fourfields.png`.

## 3. Stage 3 — the three-axis comparison   (filled after run)
