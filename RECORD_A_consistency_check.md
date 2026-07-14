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

## 2. Stage 2 — consistency verification   (filled from the full run)

## 3. Stage 3 — verdict and actions   (filled after §2)

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
