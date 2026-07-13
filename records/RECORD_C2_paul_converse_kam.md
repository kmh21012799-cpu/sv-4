# RECORD C2 — converse-KAM on Paul's critical-overlap fields (m=4,12,20,36)

**Question (Paul et al. 2022, §7):** all four critical-overlap fields have *zero*
KAM surfaces yet very different heat transport; Paul states converse-KAM (like
the Chirikov criterion) would classify all four as "very far from integrability"
— i.e. **would not distinguish them.** Nobody has tested this. We do.

**Character: verification, and there is NO answer key — only Paul's prediction.**
Outcome A (areas ≈ equal → Paul right), B (areas differ → Paul wrong, *suspect
the code*), and C (areas equal but t_c / spatial structure differ) are **all
valid results.** **No "discovery" language — observations only.**

> **★ STATUS: preparation + diagnostics done; C2-B gate PASSED; C2-C four-field
> campaign run at N=160, t_f=200. Results + interpretation below.**

---

## LIMITATIONS (read first)

1. **No answer key for the four fields themselves.** The only ground truth is
   (a) the KMM integrable island (C1) and (b) a *single-resonance* Paul field
   (planned C2-B, not yet done). Everything on the four chaotic fields is
   descriptive.
2. **Area precision is limited by cell-counting.** On the KMM answer key the
   relative error is **non-monotone in N** (+0.27% at N=90, −0.64% at N=120,
   +0.08% at N=160) — ~±0.5% aliasing from grid cells straddling the detection
   boundary (where t_c → ∞). So a *single* field's absolute area carries ~±0.5%
   method noise; only **field-to-field differences on the same grid/protocol**
   are robust, and differences below ~0.5% cannot be resolved.
3. **[UPDATED by diagnostics — the "fragmentation" was a finite-time artifact.]**
   The prep note that m=20/36 "fragment into partial bands" was WRONG: it came
   from too-short (3000-turn) orbits. Longer integration (10^5 turns) shows a
   *single* m=36 orbit spans the full [1/8, 7/8] (see Diagnosis 2 below). The
   chaos is **connected** across the domain in all four; m=36 just transports
   *slowly* (cantori). So area alone will not cheaply separate the fields — but
   for the right reason (connected chaos + KAM destroyed everywhere), which is
   what makes t_c / spatial structure the discriminating metrics.
4. **Resonance-count discrepancy with the instruction.** The instruction's P-2
   note says "m=36 has 14 resonances"; the [1/8, 7/8] separatrix condition
   requires *unit-step* n (Δψ=1/m), giving **27** resonances for m=36 (and
   3/9/15/27 for m=4/12/20/36). Our fields use the unit-step set, which
   reproduces ρ=1/8, 7/8 exactly (verified). Flagging the mismatch;
   code-suspicion default says re-check against the paper's exact mode list in
   C2-B.
5. **Euclidean metric is a modelling choice.** Paul's (ρ,θ,ζ) are orthogonal
   with |∇ρ|=1, so ξ=∂ρ ∥ ∇ρ and g=diag(1,1,1) gives the correct sign properties
   (λ(ξ)>0, λ(V)=0); it is validated on the *integrable* single resonance (where
   condition (ii)/λ is auto-satisfied by monotone rotation, so the test is
   insensitive to the metric). For the chaotic fields the metric enters condition
   (ii); a different orthogonal metric could shift a few near-threshold points.
   The core result (100% in [0.25,0.75], flat t_c) is far from that margin.
6. **The "27 vs 14 resonances" mismatch stands.** Our unit-step set reproduces
   ρ=1/8,7/8 and Chirikov S=1.000 exactly, so it is internally correct; but it
   differs from the instruction's "14". If the paper's m=36 used a different mode
   list, absolute numbers would move (the *qualitative* verdict — indistinguish-
   able core, flat t_c — is robust to this).

---

## P-1 — grid-resolution convergence (KMM example 1, answer key S_I = 0.49277)

converse-KAM non-existence area vs analytic island area, t_f = 200, 4 cores:

| N | wall-time | S(t_f) | abs err | rel err |
|---|---|---|---|---|
| 90 | 208 s | 0.49412 | +0.00135 | **+0.27%** |
| 120 | 329 s | 0.48963 | −0.00314 | **−0.64%** |
| 160 | 564 s | 0.49316 | +0.00039 | **+0.08%** |

- **N=160 reaches <0.1%** (matches the paper's 160×160 choice) — **use N=160 for
  C2.**
- The error is **non-monotone** (aliasing at the detection boundary), so read
  N=160's 0.08% as "≲0.5% method floor," not a guaranteed 0.08%. (t_f is not the
  bottleneck: the C1 curve already saturates by t_f≈130, so t_f=200 captures the
  interior; the residual is discretisation.)

## P-2 — per-configuration cost and budget (Paul fields)

Per-point cost, probed on the plain field-line integration (to ζ=t_f=200) ×2.5
variational overhead (conservative; detected points terminate earlier):

| field | # resonances | s/point (field line) | est. s/point (converse-KAM) |
|---|---|---|---|
| m=4  | 3  | 0.038 | 0.094 |
| m=12 | 9  | 0.089 | 0.222 |
| m=20 | 15 | 0.134 | 0.335 |
| m=36 | 27 | 0.218 | 0.544 |

Total budget (4 cores, N²=25600 points at N=160):

| N | m=36 alone | all four (conservative) |
|---|---|---|
| 120 | ~33 min | **~72 min** |
| 160 | ~58 min | **~128 min** |

→ **C2 at N=160 ≈ 2 h wall-time on 4 cores** (likely less; Paul points detect
early). N=120 (~72 min) is a cheaper fallback if only relative differences are
needed.

## P-sanity — Paul critical-overlap fields (the 함정-3 check, done up front)

- **Outermost separatrices at ρ = 1/8 and 7/8: verified** for all four (outer
  half-widths 0.125 / 0.042 / 0.025 / 0.014 land the outer chains' separatrices
  at 1/8 and 7/8). The *island chains* occupy the same radial range by design.
- **Chaos connectivity differs sharply** (single long field lines):
  - m=4 : one connected sea, ψ ∈ [0.130, 0.870] ≈ [1/8, 7/8].
  - m=12: one connected sea, ψ ∈ [0.122, 0.878].
  - m=20: confined to ψ ∈ [0.280, 0.878] from a mid seed (not spanning the
    lower quarter).
  - m=36: **fragmented into partially-connected bands** — lower [0.12, 0.40],
    middle [0.19, 0.59], upper [0.64, 0.88] — with regular orbits (island
    O-points / cantori) between them; **no single orbit spans [1/8, 7/8].**
- **Interpretation (kept neutral):** this is consistent with Paul's picture
  (m=36 = locally destroyed surfaces but transport-blocked by layered partial
  barriers). It means "comparable chaos volume" is a statement about *separatrix
  location*, not about *transport connectivity*. **This is exactly why C2 must
  report t_c distribution, undetected-fraction, and spatial structure — not area
  alone.**

---

---

## DIAGNOSTICS (checkpoint 2) — done before choosing N or touching the C2 body

### Decision 1: the non-monotone ±0.5% area error → **CASE B (pure discretisation)**

Separating the two error sources on the KMM answer key (analytic invariant lets
us cell-count the *true* island {Ψ<Ψ_X} directly, with no converse-KAM):

- **The wobble is 100% discretisation.** The true-island cell-count alone
  reproduces the exact non-monotone scatter — N=90:+0.27%, N=100:−1.25%,
  N=120:−0.64%, N=130:+0.68%, N=160:+0.08%, N=180:+0.72%, … — and it does **not**
  shrink with N (still ±0.7% at N=180). It is aliasing of the thin crescent
  boundary (25–27% of island cells are boundary cells), not convergence.
  `figures/diag_grid_noise.png`.
- **converse-KAM physics is exact.** From the saved t_c maps: **0 timeout
  points** (every island point detects by t_f=200) and **0 false positives**;
  the converse-KAM area equals the true-island cell-count *to the digit*. So the
  area error is entirely geometric — Diagnosis 1 (timeout) and Diagnosis 3
  (integration tolerance) are ruled out; it is Diagnosis 2 (boundary).
- **It is reducible.** 4×4 sub-cell **offset-averaging** collapses the wobble to
  **±0.1% at every N** (`diag_grid_noise.png`, orange). Cost is k² extra
  evaluations — cheap on the analytic KMM island, but 16× on the expensive Paul
  maps (so for C2 either boundary-only refinement, or a naive-N area carried with
  a ±0.3–0.5% error bar, or lean on t_c/spatial metrics).

**Resolution (per the case-B playbook):** the area wobble is a real, *bounded,
geometric* error bar (~±0.5% naive at N≥120, ~±0.1% offset-averaged). C2 area
comparisons must be stated only outside that bar; where field differences fall
inside it, report "converse-KAM does not distinguish them at this resolution"
(a valid result) and defer to the t_c distribution and spatial structure — which
are far less boundary-sensitive.

### Decision 2: m=36 "fragmentation" → **POSSIBILITY B (cantori); code validated**

- **Amplitudes are exactly critical (Diagnosis 3/4, possibility C ruled out).**
  All four fields: half-width/(Δψ/2) = 1.000, Chirikov S = 1.000, inner/outer
  separatrices at exactly ρ = 0.1250 / 0.8750.
- **The sea is connected; transport is slow (Diagnosis 1, possibility A ruled
  out).** A single m=36 orbit's cumulative ψ-span grows with time:
  10³ turns → [0.380,0.566] (w=0.19); 10⁴ → [0.123,0.655] (0.53);
  **10⁵ → [0.122,0.878] (0.76) = [1/8,7/8].** `figures/diag_m36_span.png`.
- **Multi-IC confirms one sea (Diagnosis 2).** Seeds 0.20/0.35/0.65/0.80 (20k
  turns) give overlapping ranges with **union [0.122, 0.878] ≈ [1/8, 7/8]**;
  the 0.50 seed is trapped only because it sits on the 18/36 resonance O-point.
- **The earlier "fragmentation" was a finite-time observation error** (3000-turn
  orbits looked confined) — the same finite-time trap flagged for D1. Corrected.

**Resolution:** m=36 has connected chaos across [1/8, 7/8] with strong cantori
(expected at S=1). Code is correct; no A/C bug. **NOT claimed here:** the
downstream idea that "converse-KAM can't see cantori, so it can't separate the
four" — that is a hypothesis for the C2 body / C3, recorded in
`future_questions.md`, **not** an assertion.

---

## C2-B — Paul-coordinate setup + single-resonance GATE (PASSED)

converse-KAM in Paul coordinates (ζ = time; state (ρ=ψ, θ); tangent stays
poloidal because V^ζ = 1). Ingredients: **ξ = ∂ρ** (= ∇ρ, since the (ρ,θ,ζ)
coordinates are orthogonal with |∇ρ| = 1), **Euclidean metric** g = diag(1,1,1),
**β(η,ξ) = −η^θ**, and **λ(η) = η^ρ − V^ρ g(V,η)/g(V,V)** — the exact analogues
of the KMM construction (`tools/field_paul.py`, `tools/converse_kam.py`).

**Gate.** Single-resonance (m=4, n=2, ε=0.01) is integrable (invariant
K = χ − (n/m)ψ conserved to 2.5e-10). converse-KAM run on it (N=140) detects
**exactly** the analytic island {K < K_X}: **2488 cells detected = 2488 cells
inside, 0 false positives, 0 missed** (`figures/paul_single_ckam.png`). The
converse-KAM area equals the same-grid island cell-count to the digit (0.8091);
the 1.5% vs the fine answer key (0.7973) is pure N=140 discretisation of this
wider geometry. **The Paul-coordinate converse-KAM reproduces the analytic
island — setup validated.**

## C2-C — the four critical-overlap fields (N=160, t_f=200)

t_f-convergence (N=60, t_f=200/400/800): area changes by ≤0.5% and the
undetected fraction is flat → **t_f=200 sufficient; the undetected points are
genuine non-detection, not timeout** (함정 2 discharged — undetected is *not*
read as "KAM present": raising t_f four-fold does not detect them).

| field | # res | area (%domain) | core ρ∈[0.25,0.75] | undetected % | t_c median | t_c IQR |
|---|---|---|---|---|---|---|
| m=4  | 3  | 65.9 ± ~1 | **100.0%** | 33.7 | 18.9 | [14.4, 29.8] |
| m=12 | 9  | 71.9 ± ~1 | **99.9%**  | 27.6 | 19.5 | [14.5, 30.6] |
| m=20 | 15 | 73.4 ± ~1 | **100.0%** | 26.1 | 19.4 | [14.8, 34.2] |
| m=36 | 27 | 74.0 ± ~1 | **100.0%** | 25.5 | 20.5 | [14.2, 33.6] |

(Error bar ~±1%: the single-resonance gate showed ckam==cell-count exactly, so
the area uncertainty is the geometric discretisation measured there / in the KMM
scan.) Figures: `paul_ckam_maps.png` (detection maps), `paul_tc_hist.png`
(overlaid t_c), `paul_ckam_summary.png`.

**Paul's transport ordering (Fig. 5, read qualitatively — exact numbers are
figure-only, not digitised here):** the effective parallel-diffusion volume
V_PD *decreases* with m, and the temperature flattening (heat leak) is strongest
at m=4, weakest at m=36 (m=36 insulates). i.e. transport ordering
m=4 > m=12 > m=20 > m=36.

### What converse-KAM sees (three readings, all pointing the same way)

1. **Paul's normalisation band ρ∈[0.25,0.75]: 100.0% non-existence for ALL four**
   (m=12 is 99.9%). In the transport-relevant core, converse-KAM gives the
   *identical* verdict for every field — **it cannot distinguish them.**
2. **t_c distributions are near-identical** — medians 18.9/19.5/19.4/20.5, the
   overlaid histograms lie on top of each other (`paul_tc_hist.png`). So the
   "결말 C" escape hatch (area equal but t_c differs) does **not** apply: t_c
   carries no separating signal either.
3. **Full-domain area does differ** (65.9→74.0%) — but it (a) **saturates by
   m=12** (m=12/20/36 within ~2%, near the error bar), (b) is driven by the
   *island-core / outer-edge* geometry (the detection maps show the interior
   [1/8,7/8] fully covered for all four; the difference is m=4's large regular
   island cores + coarse scalloped edge, i.e. the §4.2 island-ambiguity, not the
   chaotic sea), and (c) is **ANTI-correlated with transport** — m=4 has the
   *smallest* non-existence area yet the *largest* transport (`paul_ckam_summary.png`).

### Verdict — essentially outcome A (Paul's prediction holds)

**converse-KAM does not distinguish the four fields in any transport-relevant
way.** The core-band verdict is identical (100% for all); the t_c distributions
are identical; the only measurable difference (full-domain area) saturates by
m=12 and points *opposite* to transport, tracing island-core area rather than
the chaotic transport that V_PD measures. This is the **first direct test** of
Paul's 2022 statement that converse-KAM would not separate these fields, and it
confirms it.

**Code-suspicion check (함정 3).** The one difference is *not* the dangerous
direction (m=4 came out with the *smallest* area, not the largest), so it does
not spuriously imply "converse-KAM predicts transport." It is consistent with
(i) the validated single-resonance gate, (ii) the C1 island-ambiguity limitation
(§4.2), and (iii) exactly-critical amplitudes (Chirikov S=1.000, separatrices at
1/8,7/8, verified in P-sanity). No result here is called a "discovery" — these
are observations.

**Why (hypothesis, NOT asserted — for C3):** converse-KAM detects only the
absence of KAM surfaces; the structure that throttles m=36's transport is
*cantori* (partial barriers, confirmed in Decision 2), which are not KAM
surfaces and to which converse-KAM — and its t_c — are blind. A stickiness
diagnostic (WBA) and V_PD respond to cantori; converse-KAM does not. This is the
C3 three-axis question and is left to `future_questions.md`.

## Reproduce

```
python3 scripts/precision_c2.py sanity p2cost p1      # prep: sanity + cost + grid conv
python3 scripts/diagnose_c2.py                        # checkpoint-2 diagnostics
python3 scripts/paul_ckam_validate.py --N 140         # C2-B gate (single resonance)
for k in m4 m12 m20 m36; do
  python3 scripts/paul_ckam_campaign.py $k --N 160 --tf 200   # C2-C main
done
python3 scripts/paul_ckam_analysis.py table hist maps conv    # table + figures
python3 scripts/paul_summary_fig.py                            # summary figure
```
