# future_questions.md

Deferred questions raised by C0/C1. **Kept out of scope on purpose** — recording
them so they are not silently dropped, and so the C2/C3 work has a starting list.
None of these are claims; they are open items.

## Directly owed follow-ups (verification debt)

- **Port / cross-check the QUASR residue routine.** `tools/residue.py` is an
  independent re-implementation (QUASR not in session scope). When QUASR is
  available, cross-check R_O, R_X to 3 digits on the KMM (2,1) chain, and
  confirm the device-104183 pipeline agrees.
- **Pixel-level comparison to the source figures.** C0/C1 match Fig. 4/5/6
  *structurally* and *in trend*. Overlaying the actual published figures would
  turn "qualitative match" into a quantitative one.
- **`t_c ~ (π/2)T/√R` constant.** Measured min-`t_c` runs ~0.75–0.95 of the
  asymptotic estimate and the ratio drifts with ε (pendulum approximation
  degrading). Worth deriving the exact prefactor and the leading correction.

## The island/chaos indistinguishability (the crux for C3)

- **converse-KAM with ξ = ∇ψ cannot separate "torus of another class (island
  interior)" from "chaos"** — both are flagged non-existent (paper §4.2, and
  reproduced here: in the integrable example 1 the *island interior* is flagged,
  with no chaos present). The paper's proposed fix is a **foliation built on the
  island's elliptic field-line centre** instead of the radial ∇ψ. Implementing
  that alternate direction field would let converse-KAM distinguish island tori
  from chaotic non-existence.
- **This is exactly the axis on which converse-KAM and WBA are expected to
  disagree**: converse-KAM marks an island as "dead" (no radial torus), a
  regularity diagnostic marks it "alive" (regular motion). The disagreement set
  *is* the island set. Quantifying that overlap/disagreement is a C3 question,
  **not touched here.**

## Toward C2 (apply to Paul's critical-chaos fields)

- Run converse-KAM on `paul_m4/m12/m20/m36`. Nobody has. Expected: the
  non-existence region grows with m even as transport (heat leak) falls — the
  "broken ≠ leaky" story. Needs the Paul-field Jacobian + adapted metric
  (analogue of what was added for KMM).
- Which direction field ξ is natural for the Paul geometry (ρ,θ,ζ orthogonal)?
  ∇ψ = ∇ρ is the obvious first choice; the island-centre foliation is the
  refinement.

## Cantori hypothesis (raised by C2 diagnostics — NOT a claim)

- C2 diagnostics showed the m=36 critical-overlap field has chaos **connected**
  across [1/8, 7/8] (a single orbit spans it by ~10^5 turns) but with **slow
  transport** — strong cantori (partial barriers), as expected at Chirikov S=1.
- **Open hypothesis (still NOT asserted):** converse-KAM detects only the
  *absence of KAM surfaces*, and cantori are not KAM surfaces, so converse-KAM is
  likely blind to the very structures (cantori) that throttle m=36's transport.
- **C2 result is CONSISTENT with it** (RECORD_C2): converse-KAM gives an
  identical verdict for all four fields in Paul's core band (100% non-existence),
  identical t_c distributions, and a full-domain-area difference that runs
  *opposite* to transport — i.e. it does not track V_PD, exactly as it would if
  it were blind to the transport-controlling cantori. **This is consistency, not
  proof of the mechanism.** The mechanism claim needs C3: add V_PD and a
  stickiness diagnostic (WBA) on the *same* field and show they separate what
  converse-KAM cannot. Only then is it more than a hypothesis.

## Toward C3 (V_PD vs converse-KAM — Paul's "future publication")

- Implement the effective parallel-diffusion volume `V_PD` (Paul 2022) and
  compare, in the *same* field, with the converse-KAM non-existence volume, in
  the small-perpendicular-diffusion limit Paul names. Then add a regularity
  diagnostic (WBA) as the third axis.
- Open: does `V_PD` agree with the converse-KAM volume, the island-resolved
  version, or neither? The island/chaos distinction above is likely where the
  three metrics separate.

## Numerical / method items

- **Symmetry-line over-sampling (Thm 3.2).** θ=0 passes through every primary
  island's elliptic point, so it over-represents islands and cannot estimate
  area — confirmed here. Area must come from the 2-D grid (Thm 3.1). This is the
  same θ₀=0 trap seen before (D1/B2). Any future speed-up via symmetry must keep
  a 2-D grid for area.
- **`killends` extension** (orbits leaving the domain) — left for future work by
  the paper; not implemented.
- **Guiding-centre extension** — not implemented (future work).
- Adaptive grid refinement near separatrices would sharpen the area estimate at
  fixed cost (t_c → ∞ at the separatrix makes the boundary the slowest part).
