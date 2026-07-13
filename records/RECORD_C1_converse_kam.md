# RECORD C1 — converse-KAM 3D (MacKay 2018 / KMM 2023), validated on the integrable KMM island

**Scope (fixed in advance):** implement the converse-KAM non-existence test and
*validate it against known answers* on the KMM fields. **No new-discovery
claims, neutral outcome.** Paul-field application, V_PD, killends,
guiding-centre, and WBA comparison are all explicitly out of scope here.

Only proceeded because **C0 passed** (`RECORD_C0_fields.md`).

---

## LIMITATIONS (read first)

1. **converse-KAM with ξ = ∇ψ cannot distinguish an island interior from
   chaos.** Both are flagged "no torus transverse to ξ" (paper §4.2). Here, in
   the *integrable* example 1, the flagged region is the *island interior* (no
   chaos is present), so this is benign — but it is exactly the ambiguity that
   matters for C3 (converse-KAM marks an island "dead"; a regularity diagnostic
   marks it "alive"). The paper's fix (a foliation on the island's elliptic
   centre) is **not implemented** — deferred to `future_questions.md`.

2. **The symmetry-line "factor of 2" was only PARTIALLY reproduced.** A literal
   "check condition (i) only" reading of Theorem 3.2 does **not** halve the
   first-detection time in my implementation — on the θ=0 line the β sign-change
   (condition i) occurs at the *full* time. The factor-of-2 shows up instead in
   the **λ-condition timing**: on the θ=0 line the first `λ<0` occurs at ≈ half
   the first β-sign-change time (numbers below). I report this honestly rather
   than forcing agreement; **code-suspicion is the default** and the exact
   Theorem-3.2 detection criterion may be subtler than implemented. What *is*
   solid: the symmetry short-cut over-detects when misapplied to a 2-D grid, so
   area must come from the full Theorem-3.1 grid (shown).

3. **`t_c ~ (π/2)T/√R` is asymptotic.** Measured min-`t_c` runs ~0.75–0.95 of
   the estimate, drifting with ε (pendulum approximation). The *shape* (`∝
   1/√R`) is reproduced; the constant is not exact.

4. **Grids are 90–160, finite.** Area estimates carry O(1/N) discretisation
   error (the example-1 area overshoots S_I by ~0.3% at N=90 — a cell-size
   effect at the separatrix, where `t_c → ∞`).

5. **QUASR residue not ported** (as in C0) — the residue used in the `t_c`
   cross-check is the independent re-implementation.

---

## Method (`tools/converse_kam.py`)

Integrate the field line **and** an infinitesimal displacement vector η
(variational equation) along the KMM V-flow (φ = time), starting η(0) = ξ = ∂ψ:

```
d(ψ,ϑ)/dφ = (V^ψ, V^ϑ),      dη/dφ = J η,   J = ∂(V^ψ,V^ϑ)/∂(ψ,ϑ)  (analytic)
```

- direction field **ξ = ∂ψ** (∝ ∇ψ in the adapted diagonal metric);
- **β(η,ξ) = −η^ϑ** on the poloidal section (because V^φ = 1);
- adapted metric `ds² = dψ²/(2B₀ψ) + (2ψ/B₀)dϑ² + R₀²dφ²`;
- `λ(η) = η^ψ − V^ψ g(V,η)/g(V,V)` (the R² factor cancels; `λ(V)=0`, `λ(ξ)>0`).

**Theorem 3.1:** non-existence of an invariant torus transverse to ξ through the
orbit once, at some φ, **(i)** β changes sign **and** **(ii)** λ(η) < 0. The
whole orbit is then in the non-existence region. **Theorem 3.2** (IC on a
symmetry line): condition (i) alone.

Detection uses a terminal event on the η^ϑ zero-crossing (island orbits stop
early); a short warm-up moves η^ϑ off its exact zero start. `t_c` per grid point;
row-parallel over cores.

---

## Validation 1 — area convergence to the analytic island (Fig. 4)

Example 1 (ε=0.004), 90×90 grid in symplectic coords, `t_f = 200`:

- **`S(t_f=200) = 0.4941`  vs analytic  `S_I = 0.4928`  → ratio 1.003.**
- `figures/ck_ex1_area_convergence.png`: `S(t_f)` rises **monotonically** from
  first detection (~t_f≈18) and saturates at `S_I`. This reproduces paper Fig. 4.
- `figures/ck_ex1_map.png`: the non-existence region is the (2,1) island — two
  wide "banana" crescents at r≈0.5 (each island spans ϑ∈(−π/2,π/2)); the core
  (r<0.37) and exterior (r>0.62) are **not** flagged (radial KAM tori survive
  there); `t_c → ∞` at the separatrix (slow, bright rim near the X-points).

**This is the central C1 result: the converse-KAM non-existence area recovers the
known island area.**

## Validation 2 — first-detection time vs Greene residue (Fig. 5)

`min t_c` over the island vs `(π/2) T/√R`, `T = 4π`:

| ε | R_O | measured min t_c | (π/2)T/√R |
|------|------|------|------|
| 0.002 | 0.284 | 35.15 | 37.06 |
| 0.003 | 0.408 | 28.48 | 30.90 |
| 0.004 | 0.521 | 24.49 | 27.35 |
| 0.006 | 0.710 | 19.71 | 23.43 |
| 0.008 | 0.850 | 16.83 | 21.41 |
| 0.010 | 0.943 | 14.85 | 20.33 |

`figures/ck_ex1_tc_residue.png`. The `1/√R` law is reproduced; the measured
constant is ~0.75–0.95× the asymptotic estimate (see LIMITATION 3). This ties
the converse-KAM detection time to the independent residue computation.

## Validation 3 — symmetry short-cut vs full grid  (the D1/B2 trap)

`figures/ck_ex1_sym_vs_grid.png`, 100×100:

- **Area, Thm 3.1 (i+ii) = 0.4866; Thm 3.2 short-cut (i only) = 0.4866 —
  identical (difference +0.0%).** For the *integrable* example the tangent vector
  rotates monotonically, so the first β sign-change always occurs where η ≈ −ξ
  and hence λ<0 already holds — condition (ii) is redundant there. So the two
  criteria give the *same* non-existence region. (Over-detection by the short-cut
  is a *chaotic*-case phenomenon, not present here.)
- **Where the "factor of 2" actually lives — the λ timing.** On the θ=0 line:

  | ψ | t(β sign change) | t(first λ<0) | ratio |
  |------|------|------|------|
  | 0.10 | 31.33 | 19.52 | 0.62 |
  | 0.13 | 24.53 | 12.40 | 0.51 |
  | 0.16 | 29.19 | 12.88 | 0.44 |

  The first `λ<0` occurs at ≈ **half** the first β-sign-change time — consistent
  with the paper's statement that the symmetry-line calculation detects in half
  the time, *if* symmetry-line detection is tied to the λ-condition (quarter
  libration) rather than to the β sign-change (half libration). My literal
  "condition-(i)-only" short-cut does **not** by itself halve `t_c` (LIMITATION
  2).
- **Over-sampling (the D1/B2 trap).** The θ=0 line passes through **both**
  O-points of the (2,1) chain (and would pass through every primary chain's
  elliptic point). A 1-D scan on θ=0 therefore hits every island but sees none of
  its perpendicular extent — so **area must come from the 2-D grid**, exactly as
  the paper (and the earlier D1/B2 experience) warns.

## Validation 4 — examples 2 and 3 (both chains detected)

120×120, `t_f=200` (`figures/ck_ex2_map.png`, `figures/ck_ex3_map.png`):

- **ex 2 (2/1 + 3/2):** one broad non-existence annulus, r ≈ 0.37–0.72, covering
  *both* chains **and the chaotic layer between them** (area 0.932). Because the
  chains are bridged by chaos, converse-KAM (ξ=∇ψ) flags the whole band — it
  cannot separate the two chains or the chaos (the §4.2 limitation, in action).
- **ex 3 (2/1 + 5/4):** **two distinct** non-existence rings — the (2,1) chain at
  r ≈ 0.51 and the (5,4) chain at r ≈ 0.74 (5-fold/pentagonal structure) —
  separated by a gap of **surviving** KAM surfaces (area 0.800). Here the small
  ε leaves KAM tori between the chains, so both chains are detected *separately*.

Both examples satisfy "each chain is detected as a non-existence region."

---

## C1 pass criteria — status

| Criterion | Status |
|---|---|
| ex 1 area `S(t_f) → S_I` (Fig. 4) | **PASS** — ratio 1.003 |
| `t_c ~ (π/2)T/√R` matches residue (Fig. 5) | **PASS (shape)** — 1/√R reproduced, constant ~0.8× |
| examples 2, 3: both island chains detected | **PASS** — ex3 shows two separate rings; ex2 shows both (merged through chaos) |
| symmetry line vs 2-D grid difference recorded | **PARTIAL** — coincide for integrable ex1; factor-of-2 only via λ-timing (LIMITATION 2) |

**Verdict: PASS on the two answer-key reproductions (Fig. 4 area, Fig. 5 t_c),
PARTIAL on the symmetry-line factor-of-2 (recorded, code-suspicion noted).**

---

## Reproduce

```
python3 scripts/converse_kam_kmm.py tcres                 # Fig. 5 cross-check (fast)
python3 scripts/converse_kam_kmm.py ex1  --N 90  --tf 200 # Fig. 4 + detection map
python3 scripts/converse_kam_kmm.py sym  --N 100 --tf 200 # symmetry vs grid
python3 scripts/converse_kam_kmm.py ex23 --N 120 --tf 200 # examples 2, 3
```
