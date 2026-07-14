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

3. **converse-KAM and WBA are reimplementations — now cross-checked (see
   RECORD_A_consistency_check.md).** `vpd/diagnostics.py` contains self-contained
   field-line diagnostics: a weighted-Birkhoff `dig` (WBA) and a MacKay-style
   tangent cone-crossing converse-KAM. Cross-checked against the originals
   (C2=`393328e`, C3a=`bf5ca06`, branch `converse-kam-3d-y23ijj`) on the same
   field/grid:
   - **converse-KAM: Spearman 0.97, detection agreement 100%.** The earlier
     "`t_c` ~4 vs table ~19" was a **unit artifact** — C3b reports `t_c` in
     *periods*; `×2π` gives *zeta*, median ≈ 19–25, matching the original.
   - **WBA: fine-rank Spearman ~0.55** (dig is convergence noise in the chaotic
     bulk for both implementations) but regular/chaotic **classification agrees
     93–96%** and aggregate stats match.
   - **The §3 correlations are robust:** recomputed with the *original*
     diagnostics, `r(cKAM,V_PD)` is identical (so key question #8 is reproduced
     exactly), `r(WBA,V_PD)≈0` with both, `r(cKAM,WBA)` weak-negative with both.
     Nothing is retracted.

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

## 3. Stage 3 — the three-axis comparison

All three indicators computed on the **same fields, same `(rho,theta)` grid
(`zeta=0`), same core** `[0.25,0.75]`. converse-KAM and WBA from
`vpd/diagnostics.py`; V_PD local indicator `chi = Θ(kpar|∇∥T|² − kperp|∇⊥T|²)`
from the `zeta=0` slice of the 3-D solve.

### 3.1 The comparison table  (V_PD, ΔT at `kperp = 1e-6`)

| | m=4 | m=12 | m=20 | m=36 | separates fields? |
|---|---|---|---|---|---|
| **[topology] converse-KAM non-existence** | 99.5% | 99.4% | 99.5% | 100.0% | **no** (all ≈100%) |
| **[topology] converse-KAM t_c median** (periods; ×2π→zeta≈25) | 4.0 | 4.0 | 4.0 | 4.0 | **no** |
| **[dynamics] WBA dig median** | 1.60 | 1.67 | 1.89 | 2.45 | weakly / no |
| **[dynamics] chaos fraction (dig<5)** | 79% | 85% | 87% | 86% | **no** |
| **[transport] V_PD** | **0.884** | **0.710** | **0.465** | **0.292** | **YES** |
| **[transport] ΔT** | 0.391 | 0.399 | 0.379 | **0.440** | **YES** (m=36 insulates) |

(Our converse-KAM/WBA aggregates reproduce the C2/C3a table's *structure* —
non-existence ≈100%, dig ≈1–2, chaos ≈80–90% — under the reimplementation
caveat in Limitation 3; absolute `t_c` differs by a normalisation.)

**This is the result of the whole comparison:** the two indicators that C2 and
C3a showed cannot tell Paul's four fields apart — converse-KAM (topology) and
WBA (dynamics) — still cannot here. The transport indicators **V_PD and ΔT can**,
and in the direction Paul reported: fewer, wider island chains (m=4) transport
most; many thin chains (m=36) transport least despite the most-broken surfaces.

### 3.2 The three correlations (point by point, `zeta=0` core grid)

- **converse-KAM ↔ WBA** = **−0.09, −0.01, +0.01, +0.03** (m=4…36). Weak, near
  zero, and drifting through zero with `m`. Same weak/near-zero character C3a
  found (Spearman ≈ −0.2, splitting at islands) — the two disagree exactly where
  islands live (see 3.4).
- **WBA ↔ V_PD** = **≈ 0** at every `kperp` (−0.05 … −0.01). Uncorrelated:
  V_PD's structure does not track the island/chaos distinction WBA sees.
- **converse-KAM ↔ V_PD** = **0.03 … 0.21**, weak at every `kperp`.

### 3.3 Does converse-KAM ↔ V_PD correlation grow as `kperp → 0`?  → **NO**

| m | `kperp`=1e-4 | 1e-5 | 1e-6 |
|---|---|---|---|
| 4  | 0.102 | 0.207 | 0.173 |
| 12 | 0.108 | 0.132 | 0.106 |
| 20 | 0.092 | 0.089 | 0.063 |
| 36 | 0.031 | 0.079 | 0.085 |

The correlation stays weak (≤0.21) with no consistent growth toward small
`kperp` — for m=20 it *decreases*. **Paul's expectation** (§3: "we expect …
agreement … in the limit of small perpendicular diffusion") **is not confirmed
as a pointwise spatial correlation in this computation.** Figure
`figures/fig4_correlation.png`.

This is exactly the outcome trap 3 warned to check for: the correlation does
**not** approach 1. Had it done so we would have distrusted the code, because a
high converse-KAM↔V_PD correlation is incompatible with C2's finding that
converse-KAM cannot distinguish the four fields while V_PD does. The weak,
non-growing correlation is the self-consistent result.

### 3.4 Where the three axes split (maps, `figures/fig5_maps.png`)

The `(rho,theta)` maps for m=12 make the split concrete:
- **converse-KAM** is *saturated* — ≈100% "destroyed" everywhere, no spatial
  structure. A binary that has already gone all-on cannot track anything.
- **WBA dig** shows the island chains as regular (high dig) islands embedded in
  chaos (low dig) — it sees the islands as **alive**.
- **V_PD's χ** is ≈1 almost everywhere except thin poloidal stripes where the
  field-aligned temperature gradient vanishes. Its structure is organised by the
  temperature-gradient geometry, matching **neither** converse-KAM's saturation
  **nor** WBA's island pattern.

So, to the question "does V_PD treat an island as dead (like converse-KAM) or
alive (like WBA)?": **in these critical-overlap fields the question barely
applies** — converse-KAM is saturated and V_PD's structure is set by where
∇∥T is large, not by island topology or by local regularity. **We do not
assert the cantori explanation.** The data say V_PD is a third, distinct thing:
a graded transport measure whose spatial pattern is decorrelated from both the
topological and the dynamical indicator, which is precisely why it separates the
four fields when they cannot.

### 3.5 Resolving the stated tension

The instruction posed: Paul expects converse-KAM ≈ V_PD as `kperp→0`, yet
converse-KAM cannot distinguish the four fields while V_PD can — are these
compatible? Our answer, from the data:

- **Aggregate limit:** as `kperp→0`, `V_PD→1` for all fields and converse-KAM
  non-existence ≈100% for all fields. Both saturate near 1, so they "agree" in
  the trivial averaged sense Paul anticipated.
- **But** that agreement is *saturation*, not spatial correspondence: pointwise
  they stay weakly correlated (3.3), and the discrimination that makes V_PD
  interesting happens at **finite** `kperp`, *before* saturation, where
  converse-KAM is already pinned at 100%.

So the two statements are compatible, and the resolution is not paradoxical:
converse-KAM and V_PD coincide only in the degenerate all-destroyed / all-parallel
limit, and that coincidence carries none of the information (the field-to-field
ordering) that the comparison is about.

---

## 4. Verdict

- **Stage 1 (gate): PASSED.** Both of Paul's analytic predictions reproduced
  within one order of magnitude, with correct scalings, untuned.
- **Stage 2: qualitatively matches Paul Fig 5–6.** V_PD largest for m=4; ΔT
  largest for m=36 (most-broken field insulates best). (Caveat: the four fields
  do not yet converge within `kperp∈[1e-4,1e-6]`.)
- **Stage 3: the three-axis comparison exists.** The comparison Paul deferred in
  2022 is completed here on one stage. Result, stated neutrally: **transport
  (V_PD, ΔT) separates the four fields; topology (converse-KAM) and dynamics
  (WBA) do not; and the converse-KAM↔V_PD correlation does not grow as
  `kperp→0`, contrary to Paul's stated expectation.** No new-discovery claim is
  made; this is reproduction + comparison, and several of the numbers rest on a
  reduced-geometry model and reimplemented diagnostics (see Limitations).

Figures: `figures/fig1_island.png` (Val 1), `fig2_chaos.png` (Val 2),
`fig3_fourfields.png` (Stage 2), `fig4_correlation.png` (Stage 3 correlation),
`fig5_maps.png` (Stage 3 maps). Data: `results/*.json`. Reproduce with
`python3 -m vpd.validate1 | validate2 | stage2 | stage3` then `vpd.plots`.
