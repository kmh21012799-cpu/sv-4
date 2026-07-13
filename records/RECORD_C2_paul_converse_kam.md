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
5. **Paul Jacobian/metric not yet built.** The P-2 cost is estimated from the
   plain field-line integration × a 2.5 variational factor (conservative,
   full-t_f). The actual converse-KAM-on-Paul setup (direction field, metric,
   single-resonance validation) is C2-B.

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

## Plan for the C2 body (not yet executed)

- **C2-A:** implement/verify Paul direction field ξ (∇ρ; check ∇ψ ∥ ∂ψ in the
  orthogonal (ρ,θ,ζ) metric) and Jacobian; **validate on a single-resonance Paul
  field against the analytic W_{m,n} island area** (the Paul-coordinate answer
  key) before touching the four chaotic fields.
- **C2-B/C:** run all four at N=160, t_f=200 (with a t_f-doubling check on the
  undetected fraction — 함정 2); report **area (normalised over ρ∈[0.25,0.75],
  as Paul does), median t_c, t_c histogram (overlaid), timeout fraction, and
  detection maps**; place beside Paul's V_PD and ΔT read off Fig. 5.
- **Verdict A/B/C** stated with the precision caveat (LIMITATION 2) and the
  connectivity caveat (LIMITATION 3) explicit.

## Reproduce (preparation)

```
python3 scripts/precision_c2.py sanity p2cost   # boundary/connectivity + cost
python3 scripts/precision_c2.py p1              # KMM grid convergence
```
