# RECORD C3a — WBA on Paul's critical-overlap fields (m=4,12,20,36)

**Question:** C2 showed converse-KAM does **not** distinguish the four fields, yet
Paul's V_PD and heat transport do. Does the **Weighted Birkhoff Average (WBA)**
convergence rate `dig` distinguish them — and in which direction?

**Character: verification, no answer key, neutral outcome.** Pre-registered
hypotheses (to prevent post-hoc rationalising): **A** m=36 dig *lower* (sticky
cantori); **B** m=36 dig *higher* (fine islands, fast mixing); **C** no
difference. All three are valid results. **If A appears it is the "pretty"
result → suspect the code three times, above all the T-convergence.**

> **★ STATUS: P-1 (WBA implemented) + P-2 (gate) + cost done. Reporting and
> stopping before the C3a body (grid/T/scope pending user decision).**

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
4. **Noise floor ≈ 9–15 (rtol=1e-11).** Regular orbits saturate there, not at
   ∞; "regular" means dig ≳ 9, "chaotic" means dig ≲ 2 — a wide, clean gap (gate).

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

## Proposed C3a body (reduced, pending approval)

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

## Reproduce (prep)

```
python3 scripts/wba_gate.py --N 72 --T 300      # single-resonance gate map
```
