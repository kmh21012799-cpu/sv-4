# RECORD C0 — Two magnetic fields (KMM, Paul) + Poincaré/residue validation

**Scope of this stage (fixed in advance):** reproduce and *validate* two model
magnetic fields. This is a reproduction/verification exercise — **no new-discovery
claims, neutral outcome.** converse-KAM, V_PD, and QUASR fields are explicitly
out of scope here.

---

## LIMITATIONS (read first)

1. **QUASR residue code was NOT ported.** The instruction asked to import the
   Greene-residue routine from the QUASR pipeline (validated to 3 digits on
   device 104183). That repository is not in this session's scope
   (`kmh21012799-cpu/sv-4` only), so `tools/residue.py` is an **independent**
   re-implementation of the standard Greene construction. It is validated only
   *internally* (against the analytically known O/X points of the integrable
   example 1). A future cross-check against the QUASR routine is still owed.
   See `tools/SOURCE.txt`.

2. **"Qualitative match to Paul Fig. 6" is a visual judgement, not a metric.**
   I did not have the paper's figure pixels to overlay. The match claim rests on
   the *structural* features that must appear (connected chaotic sea spanning the
   interior, surviving KAM curves at the boundaries where the perturbation
   vanishes, primary chains at ψ = n/m, secondary chains inside the sea) — all of
   which are present — not on pixel agreement.

3. **Paul amplitudes are set by the pendulum (Chirikov) half-width**, i.e.
   `W = half the resonance spacing`. This is the leading-order island width; the
   true finite-amplitude separatrix is slightly narrower. So the fields are at
   *nominal* critical overlap, not tuned so that a measured last-KAM-surface sits
   exactly at criticality. This is faithful to the paper's stated prescription
   but is an approximation to "exactly critical".

4. **Chaotic seas are sampled, not proven.** "KAM destroyed in the interior" is
   asserted from a single field line wandering across the whole interior
   ψ-range (see below). That is strong numerical evidence, not a converse-KAM
   proof — proving non-existence is precisely the job of stage C1.

5. **Residue vs Fig. 5 is a trend comparison.** `R_O(ε)` and `R_X(ε)` reproduce
   the expected *shape* (see below); I did not read exact numbers off Fig. 5 to
   compare digit-by-digit.

---

## C0-A. KMM field (Kallinikos–MacKay–Martínez-del-Río 2023, arXiv:2304.09613)

**Implementation** (`tools/field_kmm.py`): the auxiliary vector field **V** with
`B = (B0 R0/R²) V`. Because `V^φ = 1`, φ is a time variable and the field line
flow is the 1.5-DOF system `dψ/dφ = V^ψ`, `dϑ/dφ = V^ϑ`. Standard parameters
`w1 = 1/4, w2 = 1, B0 = 1, R0 = 2, ζ = 0`, `f(ψ) = ψ − 4` for every mode.
Symplectic plotting coordinates `ỹ = √(2ψ) cosϑ`, `z̃ = √(2ψ) sinϑ`
(area = toroidal flux).

### Example 1 (integrable, single 2/1) — the answer key for C1

- **Invariant** `Ψ = −nψ − m A_φ` conserved along field lines to **6.3e-12**
  (max−min over a 50-turn orbit). Confirms integrability.
- **Poincaré section coincides with the Ψ level sets** — see
  `figures/kmm_ex1_poincare_invariant.png` (black = section, red = Ψ contours).
  A single (2,1) island chain (two islands) sits with O-points on the ϑ=0/π
  axis at r = √(2ψ_O) = **0.515**, X-points on the ϑ=π/2 axis.
- **Fixed points (ε=0.004):** ψ_O = 0.13247 (r=0.5147), ψ_X = 0.11747, from
  `dα/dφ = 0` on `sinα = 0` (α = 2ϑ − φ).
- **Island area** `S_I = ∫ Δψ dϑ` computed two independent ways (separatrix
  width integral vs. direct area count of `{Ψ < Ψ_X}`):

  | ε | S_I (width integral) | S_I (grid count) |
  |------|------|------|
  | 0.002 | 0.35030 | 0.35030 |
  | 0.003 | 0.42789 | 0.42789 |
  | 0.004 | **0.49277** | **0.49277** |
  | 0.006 | 0.60024 | 0.60023 |

  The two methods agree to ~5–6 digits, and `S_I ∝ √ε` (pendulum scaling), as
  they must. **S_I(ε=0.004) = 0.4928 is the target the C1 non-existence area
  must converge to** (paper Fig. 4).

### Example 2 (2/1 + 3/2) and Example 3 (2/1 + 5/4)

- `figures/kmm_ex2_poincare.png`: the (2,1) chain (2 islands, r≈0.51) and the
  (3,2) chain (3 islands, r≈0.65) with a chaotic layer between — as specified.
- `figures/kmm_ex3_poincare.png`: the (2,1) chain and the (5,4) chain (5 islands,
  r≈0.74).

### Greene residue cross-check (cf. Fig. 5)

The (2,1) chain is a **period-2** orbit of the Poincaré map (ϑ advances by ~π per
toroidal turn at resonance). Residues via the tangent (monodromy) map over 2
turns, `R = (2 − tr M)/4`:

- ε=0.004: **R_O = +0.521** (elliptic, 0<R<1), **R_X = −0.692** (hyperbolic, R<0)
  — correct signs.
- `figures/kmm_ex1_residue.png`: `R_O(ε)` rises linearly from 0
  (`R_O/ε ≈ 130–150` at small ε, i.e. `R_O ∝ ε`), reaches the period-doubling
  threshold `R_O → 1` near ε≈0.013; `R_X(ε)` is monotone negative. This is the
  qualitative structure of Fig. 5.
- This residue is also the C1 first-detection-time input: `t_c ~ (π/2) T/√R`,
  with `T = 4π` (the period-2 island-centre period), giving `(π/2)T = 2π² ≈
  19.74`.

---

## C0-B. Paul field (Paul–Hudson–Helander 2022, arXiv:2108.06328)

**Implementation** (`tools/field_paul.py`): field-line Hamiltonian
`χ = (ι'/2)ψ² + Σ ε_mn ψ(ψ−ψ̄) cos(mθ−nζ)`, flow `dψ/dζ = −∂χ/∂θ`,
`dθ/dζ = ∂χ/∂ψ`, normalisation `ψ̄ = ρ̄ = L_ζ = ι' = 1` so `ψ = ρ ∈ [0,1]`.
The perturbation carries the factor `ψ(ψ−1)`, which **vanishes at both
boundaries** — no KAM destruction there.

**Critical-overlap amplitudes (Sec. 7).** Resonances n/m are equally spaced
(Δψ = 1/m); setting the pendulum half-width equal to half the spacing gives the
closed form used here:

```
ε_mn = 1 / (16 · n · (m − n))        (ψ̄ = ρ̄ = ι' = 1)
```

Verified: for m=4, n=2 this gives ε = 1/64 and half-width `W = 2√(ε|ψ(ψ−1)|) =
0.125 = Δψ/2` exactly; the outermost separatrices land at ρ = 1/8 and 7/8 as the
paper specifies.

**Validation** (`figures/paul_m{4,12,36}_poincare.png`, section ζ=0, plotted in
(θ, ψ)):

- **m=4:** three primary chains at ψ = 1/4, 1/2, 3/4; a **single field line
  seeded near an inter-chain X-point wanders across the entire interior
  ψ ∈ [0.13, 0.87]** (99.9% of its returns), i.e. the stochastic layers merge
  into one connected sea — domain-spanning KAM surfaces are destroyed. KAM
  curves survive at ψ<1/8 and ψ>7/8. Secondary chains are visible in the sea.
- **m=12:** nine primary chains (ψ = 2/12 … 10/12); orbits from three different
  inter-chain seeds all span ψ ≈ [0.12, 0.88] → connected sea, KAM destroyed;
  larger primary islands than m=4, plus visible secondary chains.
- **m=36:** a dense, fine-grained chaotic sea filling ψ ≈ [0.1, 0.9] with many
  tiny embedded island chains — "visually chaotic yet fine-scale", the setup for
  Paul's central C2 result (KAM destroyed for both m=4 and m=36, yet m=36 barely
  leaks heat).

The **m=20** field is implemented (`paul_m20()`) but not plotted here (the three
above already establish the trend); it is available for C2.

---

## C0 pass criteria — status

| Criterion | Status |
|---|---|
| Example-1 Poincaré coincides with invariant Ψ contours (integrability) | **PASS** (drift 6e-12; figure) |
| Island area S_I computed analytically | **PASS** (0.4928 at ε=0.004; two methods agree to 6 digits) |
| Paul m=4/12/36 Poincaré qualitatively matches Fig. 6 | **PASS (visual)** — connected sea + surviving boundary KAM + primary/secondary chains |
| Residue compared to Fig. 5 | **PASS (trend)** — correct signs, R_O∝ε, period-doubling near ε≈0.013 |

**Verdict: C0 PASS** (with the QUASR-residue-port caveat and "visual"/"trend"
qualifiers recorded above under LIMITATIONS).

---

## Reproduce

```
pip install numpy scipy matplotlib
python3 scripts/poincare_kmm.py                       # KMM ex 1/2/3 figures
python3 scripts/residue_kmm.py                        # residue vs eps figure
python3 scripts/poincare_paul.py m4  --nseed 45 --ncross 200 --nchaos 30000
python3 scripts/poincare_paul.py m12 --nseed 45 --ncross 180 --nchaos 20000
python3 scripts/poincare_paul.py m36 --nseed 40 --ncross 180 --nchaos 7000
```
