# future_questions.md — items isolated out of C3b (updated after the v2 rerun)

These are questions raised but deliberately **not** answered, to keep the record
to reproduction + comparison. None is a claim.

**RESOLVED by the v2 rerun (RECORD_C3b_v2_vpd.md):**
- *Field model.* C3b originally used a constant-amplitude field missing Paul's
  `psi(psi−psibar)` envelope. Redone on Paul's actual field (original code); the
  conclusions held and sharpened.

**Still open:**

1. **Toroidal geometry.** Still the Euclidean (ρ,θ,ζ) slab — the same choice
   C2/C3a made, and Paul's constant-Jacobian metric makes the volume element
   uniform (R-weighting changes V_PD <0.05%). A full toroidal metric for the
   *operator* (not just the volume) remains untested.

2. **κ⊥ → 0 convergence.** On Paul's field the four V_PD are still separated at
   κ⊥=1e-6 (0.897→0.423); they converge to 1 only at smaller κ⊥. Needs
   κ⊥ ≤ 1e-7, where AMG already struggles (m=36 1e-6 hit the 2000-iter cap at
   residual 1e-6). A field-aligned preconditioner would be the enabling step.

3. **A graded converse-KAM.** converse-KAM is *fully saturated* (100%
   non-existence) on all four fields, so its pointwise correlation with V_PD is
   undefined (binary) or weak (graded t_c). A graded topological observable
   (largest surviving invariant set; distance-to-a-torus) might correlate where
   the saturated binary cannot — the cleanest way to actually test Paul's
   small-κ⊥ agreement expectation.

4. **Cantori.** Still not tested; the data neither confirm nor refute. A
   turnstile-flux calculation on these fields is the clean route.

5. **WBA T=5000 vs T=1000.** The v2 WBA maps used n_periods=1000 (dig median
   1.0–1.7); C3a's table used T=5000. The classification (chaos fraction) is
   robust to this, but the m=36 dig median (1.72 vs C3a 0.99) reflects the
   higher noise floor at fewer periods.
