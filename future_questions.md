# future_questions.md — items isolated out of C3b

These are questions raised but deliberately **not** answered in C3b, to keep the
record to reproduction + comparison. None is a claim.

1. **Toroidal geometry.** C3b uses a reduced slab. Does the ΔT ordering
   (m=36 insulates most) survive in a real toroidal annulus with `|B|` variation
   and a genuine `B·∇` operator? Paul's absolute numbers live there.

2. **The `kperp→0` convergence of the four fields.** Within `[1e-4,1e-6]` the
   four V_PD do not converge (spread 0.25→0.60→0.59). At what `kperp` do they
   actually merge toward 1, and does the *ordering* persist all the way down?
   Needs `kperp ≤ 1e-7` (anisotropy ≥1e7; AMG iteration counts already ~1800 at
   1e-6).

3. **Why is V_PD's χ organised by poloidal angle, not by island position?**
   The maps (fig5b) show χ=0 in thin θ-stripes where ∇∥T vanishes, not on the
   island chains. Is this a ζ=0-slice artifact, or a genuine geometric feature
   of the isotherm/field-line angle? A full 3-D χ-structure study would tell.

4. **Cantori.** We explicitly did **not** test whether the residual transport
   barriers in the high-m fields are cantori. The instruction warned against
   asserting this pre-data; the data here neither confirm nor refute it. A
   flux-through-cantori calculation (à la MacKay–Meiss–Percival turnstile) on
   these fields is the clean way to ask.

5. **A graded converse-KAM.** converse-KAM saturates (≈100% for all four fields),
   which is *why* its pointwise correlation with V_PD is weak. A graded
   converse-KAM observable (e.g. the size of the largest surviving invariant
   set, or a "distance to a torus") might correlate with V_PD where the binary
   cannot. Would that recover Paul's expected agreement?

6. **t_c normalisation.** Our converse-KAM `t_c` ≈ 4 vs the C2 table's ≈19.
   Reconciling the two implementations' detection thresholds would let the
   absolute `t_c` numbers be compared, not just their (flat) ordering.
