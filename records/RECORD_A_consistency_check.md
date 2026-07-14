# RECORD A — consistency check of C3b's reimplemented converse-KAM / WBA

**Character: consistency verification. No new physics claim. Neutral ending.**

C3b's three-axis comparison computed its point-to-point correlations with
*reimplementations* of converse-KAM and WBA, not the original C2/C3a code
(C3b honestly flagged this). This record cross-checks the reimplementations
against the originals and states how much the correlation results can be
trusted.

What is at stake (depends on the reimplementations):
`r(cKAM↔V_PD)`, `r(WBA↔V_PD)`, `r(cKAM↔WBA)`, and key question #8
("the cKAM↔V_PD correlation does not grow as κ⊥→0").
What is *not* at stake (each measured with its own original code): that
converse-KAM (C2) and WBA (C3a) cannot distinguish the four fields, and that
V_PD (C3b, validated) can.

---

## 0. LIMITATIONS (read first)

1. **The comparison holds the field fixed and varies only the algorithm.** We
   run the *original* converse-KAM / WBA (`consistency/orig/`, pinned to C3a tip
   `bf5ca06`, which also contains C2 `393328e`) on the *same C3b field* and the
   *same grid points* used by the C3b diagnostics, via a thin interface adapter
   (`consistency/adapter.py`). This isolates diagnostic-algorithm fidelity. It
   does **not** re-derive the correlations on the original *field* (see #3).

2. **Original grid-level data is gone.** C2/C3a `.gitignore`d their
   `data/*.npz`, so the original per-point `t_c`/`dig` arrays did not survive the
   container. We therefore regenerate them by running the original *code* (which
   did survive) on the same points — an equivalent and cleaner comparison.

3. **The C3b field is a *different, simpler* model than C2/C3a's.** C2/C3a use
   Paul's Hamiltonian with the `psi(psi−psibar)` envelope and a `dθ/dζ`
   perturbation; C3b uses a constant-amplitude reduced model with `iota=0.5+ρ`.
   Both use the **same Euclidean (1,1,1) slab metric and uniform (ρ,θ,ζ)
   volume** — the original C2 record documents that identical choice. So the
   geometry is consistent between C2/C3a and C3b; the *field spectrum/envelope*
   is not. This check validates the diagnostic algorithms on the C3b field; the
   field-model difference is logged as its own limitation, not silently folded in.

4. **No tuning, no new physics.** The only computation added is running the two
   original diagnostics on the C3b field (a re-use of existing code, reported
   here). No parameters were adjusted to improve agreement.

---

## 1. Stage 1 — why the originals were not used in C3b

Diagnosis (all verified against the repo, not assumed):

- **The original code was in a different branch, never fetched into the C3b
  session.** C2 (`393328e`) and C3a (`bf5ca06`) live on
  `origin/claude/converse-kam-3d-y23ijj`. The C3b session started from an empty
  `main` (`baf084f`, a stub README) and its scope listed only
  `kmh21012799-cpu/sv-4`; the converse-kam branch was not checked out, so at C3b
  time there was no `tools/converse_kam.py` or `tools/wba.py` to import. This is
  the proximate cause: **(c) code not present in the working tree.**
- **Grid data was also unavailable.** Even had the branch been known, C2/C3a
  `.gitignore` excluded `data/*.npz`, so no saved `t_c`/`dig` grids existed to
  load — **(a) data files absent.**
- **The field/grid differ too.** C3b's field model and `iota` profile differ
  from C2/C3a's (Limitation 3) — **(b) partial**, but secondary to (c).

So C3b reimplemented from the same source papers, exactly as C2/C3a themselves
had done (their `SOURCE.txt`: "written fresh from the two source papers"; QUASR
not in scope). The reimplementation was the honest option available; this record
now checks it was a faithful one.

## 2. Stage 2 — consistency verification

Original vs C3b diagnostics, **same C3b field, same 480 grid points**
(`rho∈[0.25,0.75]×24`, `theta∈[0,2π/m)×20`), `t_f=200`, WBA `n_periods=300`.

### 2.1 converse-KAM (original MacKay/KMM Thm 3.1  vs  C3b tangent test)

| quantity | m=12 | m=4 |
|---|---|---|
| detection agreement (detected/undetected) | **100.0%** | **100.0%** |
| Spearman `t_c` (detected-by-both) | **0.978** | **0.975** |
| Pearson `t_c` | 0.978 | 0.995 |
| `t_c` median, original (zeta) | 22.6 | 19.2 |
| `t_c` median, C3b (zeta) | 25.1 | 25.1 |
| scale ratio (C3b/orig), median | 1.11 | 1.10 |

→ **CASE 1 (fully consistent).** Ranks preserved (Spearman ≥ 0.97), detection
identical, absolute `t_c` matches within ~10–30%. The small ~1.1 systematic and
the C3b median pinned at `4·2π=25.1` come from C3b checking the tangent sign
once per period (period-granular `t_c`) vs the original's event-based timing.
Figure `figures/figA_consistency.png` (top).

### 2.2 The "t_c differs" flag was a **unit artifact**

C3b reported `t_c` in **periods** (median ≈ 4); the original and the C2 table
report it in **zeta** (median ≈ 19). `4 periods × 2π = 25.1 ≈ 19–25`. Once units
match, they agree — there is no physics discrepancy. (C3b's LIMITATIONS is
corrected accordingly.)

### 2.3 WBA (original continuous-integral dig  vs  C3b stroboscopic dig)

| quantity | m=12 | m=4 |
|---|---|---|
| Spearman `dig` (all points) | 0.574 | 0.514 |
| Pearson `dig` | 0.875 | 0.923 |
| **regular/chaotic classification agreement** (dig≷5) | **92.9%** | **96.3%** |
| chaos fraction (dig<5), original | 0.756 | 0.758 |
| chaos fraction (dig<5), C3b | 0.823 | 0.788 |
| dig median, original | 1.27 | 1.11 |
| dig median, C3b | 1.83 | 1.52 |

→ Fine-rank **Spearman 0.51–0.57 is below the 0.7 gate**, but this is *intrinsic*,
not a fidelity defect: `dig` in the chaotic bulk is a convergence-noise value
with no stable rank in **either** implementation (the scatter, `figA` bottom,
shows two clean clusters — both-chaotic and both-regular — with a scrambled
low-dig blob). The quantity that matters, the regular/chaotic **classification,
agrees 93–96%**, and the aggregate stats match. So the disagreement is confined
to the meaningless fine-ordering inside the chaotic sea.

### 2.4 Slab vs torus volume element (secondary)

Applying a large-aspect-ratio torus weight `dV = (1+0.3 cosθ) dρdθdζ` to V_PD:
`rel change = −0.01%` (m=12), `+0.03%` (m=4). **Negligible.** Consistent with
Paul's model having a constant metric Jacobian (uniform volume) and with the
original C2 using the same uniform `dρ∧dθ` weighting. The slab volume element is
not an approximation here.

## 3. Stage 3 — verdict and the decisive robustness test

The fidelity Spearman (0.97 for cKAM, 0.55 for WBA) bounds the *diagnostics*.
But what is actually at stake is the **correlations**. So we recomputed C3b's
three correlations **using the original diagnostics** on the same field/grid
(the field-line diagnostics are κ⊥-independent, so reused across κ⊥ while only
V_PD's χ is re-solved). `results/consistency_correlations.json`:

| correlation | field | using **original** | using **C3b** |
|---|---|---|---|
| r(cKAM, V_PD) @ κ⊥=1e-4/1e-5/1e-6 | m=12 | 0.151 / 0.088 / 0.088 | **identical** |
| r(cKAM, V_PD) @ κ⊥=1e-4/1e-5/1e-6 | m=4 | 0.089 / 0.194 / 0.194 | **identical** |
| r(WBA, V_PD) (weak/≈0) | m=12 | −0.05 … −0.09 | −0.01 … −0.02 |
| r(WBA, V_PD) (weak/≈0) | m=4 | −0.04 … −0.05 | −0.09 … −0.10 |
| r(cKAM, WBA) | m=12 | −0.163 | −0.099 |
| r(cKAM, WBA) | m=4 | −0.056 | −0.102 |

- **r(cKAM, V_PD) is identical** with the original (detection agrees 100%), so
  **key question #8** ("the correlation does not grow as κ⊥→0") is reproduced
  *exactly* by the original converse-KAM. **Fully robust.**
- **r(WBA, V_PD) ≈ 0** with both — the null survives the swap despite WBA's
  0.55 fine-rank fidelity, because a null correlation is preserved under
  noisy relabelling. **Robust.**
- **r(cKAM, WBA)** is weak-negative with both; the *original* gives −0.16 for
  m=12, even closer to C3a's reported −0.2 than the C3b reimplementation.

### Verdict

- **converse-KAM: CASE 1.** Consistent (Spearman 0.97, 100% detection agreement).
- **WBA: CASE 2** in effect. Fine-rank Spearman (0.55) is below the strict gate,
  but (i) the regular/chaotic classification agrees 93–96%, (ii) the aggregate
  stats match, and (iii) **every correlation that uses WBA is a weak/null value
  that both implementations reproduce.** No C3b conclusion depends on WBA's
  fine ranking.
- **The correlation results and key question #8 are trustworthy.** #8 rests on
  converse-KAM, which is fully consistent; the WBA-dependent correlations are
  nulls that survive substituting the original code. **Nothing is retracted.**

### Actions taken
- `consistency/SOURCE.txt`: reimplementation + cross-check result recorded.
- `RECORD_C3b_vpd.md` LIMITATIONS: corrected the `t_c` unit note; added the
  cross-check outcome (correlations robust, converse-KAM Spearman 0.97).
- No correlation was recomputed for the final C3b tables (they are unchanged);
  the robustness is documented rather than forcing a numeric edit that the
  cross-check shows would be identical.

## 4. Secondary — slab vs torus volume element

Paul's domain is a toroidal annulus, but his field-line-Hamiltonian model has
**no major-radius dependence**, and its metric
`ds² = dψ²/(2B0ψ) + (2ψ/B0)dθ² + R0²dφ²` has a **constant Jacobian**
`√g = R0/B0`. So the physical volume element is uniform in `(ψ,θ,ζ)` and the
`R`-weighting `dV ∝ R dρdθdζ` of a genuine torus does **not** apply — uniform-
volume V_PD is exact for the Paul model, not a slab approximation. The original
C2 code computes its non-existence area with the same uniform `dρ∧dθ` weighting,
and the C2 record explicitly logs the Euclidean-metric choice. A numerical
`R`-weighting sensitivity check is reported in §2.4.
