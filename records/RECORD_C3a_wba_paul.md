# RECORD C3a — WBA on Paul's critical-overlap fields (m=4,12,20,36)

**Question:** C2 showed converse-KAM does **not** distinguish the four fields, yet
Paul's V_PD and heat transport do. Does the **Weighted Birkhoff Average (WBA)**
convergence rate `dig` distinguish them — and in which direction?

**Character: verification, no answer key, neutral outcome.** Pre-registered
hypotheses (to prevent post-hoc rationalising): **A** m=36 dig *lower* (sticky
cantori); **B** m=36 dig *higher* (fine islands, fast mixing); **C** no
difference. All three are valid results. **If A appears it is the "pretty"
result → suspect the code three times, above all the T-convergence.**

> **★ STATUS: P-1 + gate + cost + IC-vectorised RK4 + Stage 1 (T-convergence) +
> Stage 2 focused (t_c↔dig correlation) all done. Verdict below. C3a complete;
> V_PD (C3b) out of scope.**

---

## LIMITATIONS (read first)

1. **`dig` is a convergence rate, not chaos itself.** Low dig means the weighted
   average has not settled at this T — from genuine chaos *or* from a
   finite-T/long-period orbit. This is why the **T-convergence check
   (T=1000/2000/5000) is mandatory** and gates every field-to-field claim (함정
   3): a difference that shrinks with T is a finite-time artifact, not structure.
2. **Full N=160 grid is INFEASIBLE.** WBA integrates 2·T periods per point
   (T=1000 → 2000 toroidal turns) vs converse-KAM's ~32. Measured per-point cost
   (T=1000): m=4 ~6–18 s, **m=36 ~10 s (regular) to ~60–100 s (chaotic)**. A full
   25600-point m=36 grid is ~110 core-hours. The C3a body must use **sampling +
   low-N maps** (proposed below), so absolute area-style numbers are replaced by
   *distributions over a sample*.
3. **Observable choice matters; we use h=ψ (=ρ).** Position (x=ρcosθ, y=ρsinθ)
   gave *worse* separation (chaotic dig_pos ~4–6 vs dig_ψ ~0.5–2) because the m
   islands sit symmetrically so ⟨x⟩,⟨y⟩≈0 and the relative digit is unusable.
   The bounded-coordinate worry (함정 1, from D_FL) does **not** bite here: WBA
   measures average convergence, not spread, and ψ∈[0,1] gives clean high dig for
   regular orbits (see gate).
4. **Noise floor ≈ 9–15 (adaptive, rtol=1e-11).** Regular orbits saturate there,
   not at ∞; "regular" means dig ≳ 9, "chaotic" ≲ 2 — a wide gap (gate).
5. **Stage 1 uses fixed-step RK4 (M=96), not adaptive.** It is ~200–300× faster
   (vectorised over ICs) and validated against adaptive on the core (median
   1.05 vs 1.13, corr 0.96). But it **under-resolves the regular-orbit tail** of
   high-m fields (m=36 island 12.8 vs adaptive 13.9; near-axis/boundary regular
   orbits worse). The *chaotic bulk* — which sets the core median — is captured
   correctly, so the Stage-1 verdict (indistinguishable core distributions) is
   robust; a small residual bias could only *depress* the high-dig tail of m=36,
   which would if anything favour the (rejected) hypothesis A, so the null result
   is conservative.

---

## P-1 — WBA implementation (`tools/wba.py`)

Bump weight `g(s)=C exp(−1/(s(1−s)))`, C=142.2503757771, verified `∫₀¹g = 1.0000`.
WBA computed as an extra ODE alongside the field line (ζ = time):
`dW/dζ = g(ζ/τ) h`, `WB = W(τ)/τ`, τ = 2π·(n_periods). Convergence judge over two
consecutive windows [0,τ],[τ,2τ]: `absdig=−log₁₀|A−B|`, `reldig` (relative),
`dig=max`. Observables computed: **dig_ψ** (h=ρ) and **dig_pos** (min over
h=ρcosθ, ρsinθ).

## P-2 — GATE (regular → high dig, chaotic → low dig)

**Point probes (T=1000):**

| field / point | region | dig_ψ | dig_pos |
|---|---|---|---|
| single-res, island ρ=0.45 | regular (libration) | 14.7 | 13.5 |
| single-res, circulating ρ=0.2 | regular | 9.3 | 14.7 |
| single-res, outer ρ=0.85 | regular | 9.1 | 13.8 |
| **m=4**, island O ρ=0.5 | regular | **13.9** | 13.8 |
| **m=4**, chaotic sea ρ=0.375 | chaos | **2.1** | 3.8 |
| **m=4**, chaotic sea ρ=0.44 | chaos | **0.5** | 3.6 |
| **m=36**, island ρ=0.5 | regular | **13.9** | 13.4 |
| **m=36**, chaotic sea ρ=0.35 | chaos | **1.5** | 6.0 |

→ **dig_ψ cleanly separates regular (≳9) from chaotic (≲2)** on the *real*
fields, and gives uniformly high dig on the integrable single resonance.
dig_pos separates worse (chaotic 3.6–6.0). **h=ψ is the observable of record.**
Gate map: `figures/wba_gate.png` (single-resonance dig_ψ). **GATE PASSED.**

## Cost (measured, T=1000, 1 core)

| field | # modes | regular pt | chaotic pt |
|---|---|---|---|
| m=4  | 3  | ~6 s  | ~18 s |
| m=36 | 27 | ~10 s | ~60–100 s |

Chaotic orbits dominate (~75% of the domain per C2) and are the expensive ones
(tiny steps over 2000 turns). **Full N=160 four-field campaign ≈ 100+ core-hours
— not run.**

## Speedup — IC-vectorised fixed-step RK4 (cost check ①)

WBA needs uniform sampling anyway, so a fixed step is natural and lets many
orbits integrate at once (one vectorised RHS per step). `wba_dig_rk4`
(`tools/wba.py`): **~80 ms/IC (m=4) to ~210 ms/IC (m=36)** vs 15–60 s/IC adaptive
— **≈200–300×**. Step count M (steps/period): M=16 fine for islands/chaos but
under-resolves (i) high-m fields' high-frequency forcing (m=36 needs M≳70) and
(ii) circulating orbits' long-term phase. **M=96 validated against adaptive on 30
random core m=36 ICs: median 1.05 vs 1.13, correlation 0.96** (chaotic bulk
matches; the regular-tail/island points are under-resolved by ≲a few dig — a
documented caveat, LIMITATION 5). Same M=96 for all four fields → unbiased
comparison. Early termination NOT used (dig needs both windows).

## Stage 1 — T-convergence of the core dig distribution (the decisive check)

100 **random** ICs (not grid/line) in ρ∈[0.25,0.75], same ICs across all fields
and T; RK4 M=96, h=ψ. **Median core dig vs T:**

| field | T=500 | T=1000 | T=2000 | T=5000 |
|---|---|---|---|---|
| m=4  | 1.05 | 1.19 | 1.21 | 1.46 |
| m=12 | 0.93 | 0.93 | 0.83 | 0.97 |
| m=20 | 1.21 | 1.12 | 1.01 | 0.93 |
| m=36 | 1.42 | 1.24 | 1.14 | 0.99 |

(frac(dig<5) ≈ 87–92% for every field at every T.) Figures:
`wba_stage1_Tconv.png` (median ± IQR band vs T), `wba_stage1_hist.png`
(overlaid distributions at T=5000).

### Verdict — hypothesis C (no distinction); the finite-time trap caught

- **The four core dig distributions are indistinguishable.** All sit at median
  ≈1 with ~90% of points chaotic (dig<5); the IQR bands overlap completely at
  every T (`wba_stage1_Tconv.png`).
- **No T-stable ordering.** The apparent T=500 spread (m36 highest at 1.42, m12
  lowest at 0.93) does **not** persist: as T grows **m36 falls 1.42→0.99** into
  the pack while **m4 rises 1.05→1.46**; the curves cross. A single-T snapshot
  (e.g. T=1000, m36=1.24 highest) would have mis-suggested a difference — the
  **T-convergence check prevented that**, exactly its purpose (the D1/finite-time
  trap, third encounter).
- **Hypothesis A (m36 stickier → lower dig) is REJECTED.** m36 is not the lowest;
  its short-T *elevation* is a finite-time transient that washes out. (So the
  "pretty result" did not appear — and the residual m4>rest at large T tracks
  m=4's larger island cores, i.e. is ANTI-correlated with transport, the same
  sign as C2's area — not a transport signal.)
- **WBA (core, h=ψ) does not distinguish the four fields — same as
  converse-KAM.** Both a KAM-existence indicator (C2) and an orbit-regularity
  indicator (this) fail in the core where V_PD succeeds.

**Per the plan's decision rule** (differences shrink/reshuffle with T ⇒
finite-time artifact ⇒ hypothesis A rejected ⇒ no full grid needed): Stage 1
answers the minimal proposition. Stage 2 (core reduced-grid spatial maps +
t_c↔dig correlation) remains optional and independently interesting ("what do
the two axes say about each other"); deferred to the user.

## Stage 2 (focused) — converse-KAM t_c vs WBA dig at the same points

400 random core ICs/field; C2 t_c (nearest grid cell) vs C3a dig (RK4 M=96,
T=1000). `figures/wba_tc_dig_corr.png`.

| field | core pts detected by converse-KAM | Spearman(t_c, dig) |
|---|---|---|
| m=4  | 399/400 | −0.20 |
| m=12 | 399/400 | −0.17 |
| m=20 | 399/400 | −0.21 |
| m=36 | 400/400 | −0.14 |

- **The two axes are only weakly related (|ρ|≈0.14–0.21).** Not tautological
  (|ρ|≪1), not independent (≠0) — they measure *different* things (cf. the B4
  dig↔FTLE ≈ −0.53; weaker here).
- **The scatter is L-shaped, and the disagreement is exactly the islands.** The
  chaotic bulk sits at low t_c + low dig (both axes say "chaotic"). But a
  vertical spray at low t_c carries **high dig (up to ~14)**: these are island
  interiors — **converse-KAM detects them fast ("no radial torus → dead") while
  WBA calls them regular ("alive").** This is the §4.2 converse-KAM limitation
  made quantitative, and it confirms the C1 prediction that the two indicators
  split *precisely on islands*. It is also why the correlation is weak/negative:
  agreement on the sea, anti-agreement on islands.
- converse-KAM detects ≈100% of core points (as in C2); the WBA-regular island
  points are a minority of the core but the sole locus of disagreement.

## Overall C3a verdict

**WBA does not distinguish the four Paul fields (core, h=ψ) — the same null as
converse-KAM.** Both a KAM-existence axis and an orbit-regularity axis give the
same reading for all four in the transport-relevant core, while V_PD/heat
transport separates them cleanly (Paul). The two axes are themselves only weakly
correlated and split on islands, so they are genuinely different lenses — but
*neither* is a transport proxy here. Hypothesis A (m=36 stickier → lower dig) is
**rejected** (finite-time artifact, caught by the T-scan). No discovery claim —
a clean null observation. **What is left for C3b:** implement V_PD on the same
fields and show it separates what these two axes cannot (the actual Paul
comparison); the cantori-blindness *mechanism* remains a hypothesis until then
(`future_questions.md`).

## (superseded) earlier proposed C3a body

- **dig distributions** by *random sampling* (~800–1500 pts/field, weighted to the
  core ρ∈[0.25,0.75]) at T=1000 — the primary deliverable (histograms, median,
  IQR), since "look at the distribution, not a binary chaos fraction" is the whole
  point.
- **Low-N spatial maps** (N≈48) for visualization.
- **T-convergence** (T=1000/2000/5000) on a small core sample (~200–300 pts/field)
  — the decisive finite-time check.
- **t_c↔dig correlation** at the sampled points (using the stored C2 t_c map).
- Estimated budget at T=1000 base: **~10–25 core-hours** depending on sample
  sizes; can be trimmed (T=500 base halves it; separation still clean).

## Reproduce

```
python3 scripts/wba_gate2.py                 # gate: m=4 dig map (real chaos)
python3 scripts/wba_stage1.py                # Stage 1 T-convergence (core)
python3 scripts/wba_stage2_corr.py           # Stage 2: t_c (C2) vs dig
```
