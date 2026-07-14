# converse-KAM 3D and V_PD — Comparing Topological, Dynamical, and Transport Indicators

**Part 3 of a summer toolkit.**
**Full narrative: https://github.com/kmh21012799-cpu/sv1_comp**

## The open item

Paul, Hudson & Helander (2022, *J. Plasma Phys.* **88**, 905880107) introduced
V_PD, the effective volume of parallel diffusion, and closed their paper with:

> "our metric appears to be related to the converse KAM approach... We expect
> that the effective volume of parallel diffusion might agree with such a
> calculation in the limit of small perpendicular diffusion. **We reserve such a
> comparison for future publication.**"

That comparison has not appeared. This repository carries it out.

## The test bed

Paul's model is a field-line Hamiltonian on a periodic slab (ρ, θ, ζ), ρ ∈ [0,1]:

```
χ(ψ, θ, ζ) = ι' ψ²/2  +  Σ_{m,n} ε_{m,n} · ψ(ψ − ψ̄) · cos(mθ − nζ),   ψ = ρ.
```

The `ψ(ψ − ψ̄)` envelope makes the perturbation vanish at both boundaries, so
that there is no parallel flux into the volume. The amplitudes are set to
critical island overlap, `ε_{m,n} = 1/(16 n(m − n))`, giving a Chirikov
parameter S = 1.000 and outermost separatrices at ρ = 1/8 and 7/8. Four fields
are studied, distinguished by their poloidal mode number m:

| Field | resonances n/m | number of chains |
|---|---|---|
| m = 4  | 1/4, 2/4, 3/4      | 3  |
| m = 12 | 2/12 … 10/12      | 9  |
| m = 20 | 3/20 … 17/20      | 15 |
| m = 36 | 5/36 … 31/36      | 27 |

All four fields have **no surviving KAM surfaces**: a single field line seeded in
the chaotic sea wanders across the whole interval [1/8, 7/8]. Their heat
transport nonetheless differs substantially.

## Results

Same fields, same grid, same core (ρ ∈ [0.25, 0.75]), each axis computed with its
own original code (converse-KAM from C2, WBA from C3a, V_PD from C3b) on Paul's
field. Perpendicular diffusivity κ⊥ = 10⁻⁶.

| | m = 4 | m = 12 | m = 20 | m = 36 | Distinguishes? |
|---|---|---|---|---|---|
| **[Topology] converse-KAM non-existence** | 99.7% | 100% | 100% | 100% | **No** |
| **[Topology] median t_c** (ζ) | 18.5 | 20.7 | 19.5 | 20.0 | **No** |
| **[Dynamics] median WBA dig** | 1.16 | 1.00 | 1.10 | 1.72 | **No** |
| **[Dynamics] chaotic fraction** (dig < 5) | 0.87 | 0.88 | 0.88 | 0.77 | **No** |
| **[Transport] V_PD** | **0.897** | 0.834 | 0.697 | **0.423** | **Yes** |
| **[Transport] ΔT** | **0.039** | 0.060 | 0.081 | **0.113** | **Yes** |

Commits: C2 = `393328e`, C3a = `bf5ca06`, C3b v2 = `8a52ca1`.

### Reading

**Neither topology nor dynamics separates the four fields.** converse-KAM
saturates at ~100% non-existence in the core; WBA gives median dig ≈ 1 with
overlapping distributions and a chaotic fraction near 0.8–0.9 for all four.

**Transport does.** V_PD falls monotonically from 0.897 (m = 4) to 0.423
(m = 36); ΔT rises monotonically from 0.039 (m = 4) to 0.113 (m = 36).

m = 36 has the most resonances and is the most thoroughly destroyed field. It is
also the best insulator.

## Paul's expectation does not hold

Paul et al. make two statements. One holds; the other does not.

| | Statement | Outcome |
|---|---|---|
| 1 | converse-KAM will not distinguish these fields | **Confirmed** |
| 2 | converse-KAM and V_PD will agree as κ⊥ → 0 | **Not observed** |

The pointwise spatial correlation r(−t_c, V_PD) as κ⊥ decreases (10⁻⁴ / 10⁻⁵ /
10⁻⁶):

```
m = 4:    0.116  →  0.022  →  0.043
m = 12:  −0.036  → −0.014  →  0.017
m = 20:  −0.036  → −0.041  → −0.040
m = 36:   0.002  → −0.001  →  0.007
```

The correlation does not grow toward small κ⊥; it stays within |r| ≲ 0.12 for all
four fields. The WBA–V_PD correlation is likewise ≈ 0.

**The structural reason.** converse-KAM is fully saturated in the core: it reports
100% non-existence for every field, so its spatial variance is zero. There is
nothing for V_PD to be correlated with. The two statements are in fact
incompatible — if V_PD distinguishes the fields and converse-KAM agreed with
V_PD, then converse-KAM would distinguish them too, which statement 1 denies. The
apparent "coincidence" in the small-κ⊥ limit is only that both quantities
saturate: converse-KAM → 100% and V_PD → 1 for every field.

## What this means

Topological and dynamical indicators see only the magnetic field. V_PD solves a
temperature field, and the temperature field carries the plasma state through κ⊥.
Two fields with the same magnetic structure but different κ⊥ give different
transport, and the field-only indicators cannot see that.

This is the same point made by Vlad et al. (2002): the same magnetic field can
give diffusion coefficients differing by orders of magnitude depending on the
plasma state. The connection length is only the numerator.

## Validation

Each axis was validated against a known answer before use.

**converse-KAM (C1).** On an integrable KMM field it recovers the analytic island
area to 0.3% (0.4941 vs 0.4928 at a 90×90 grid), and the first-detection time
follows t_c ∝ 1/√R (measured constant 0.75–0.95× the asymptotic estimate),
cross-checked against an independent Greene-residue computation
(R_O = +0.521, R_X = −0.692). The KMM invariant is conserved along field lines to
6.3×10⁻¹².

**WBA (C3a).** The gate cleanly separates regular from chaotic motion: island
interiors give dig ≈ 13.9 (14.7 on the integrable single resonance), the chaotic
sea gives dig ≈ 0.5–2.1 (1.5 for m = 36).

**V_PD (C3b).** Reproduces Paul's analytic predictions: the island transition at
ε_crit = √κ⊥/2 (the measured transition sits within 1.6–3.4× of this in the
effective local amplitude), V_PD ∝ √ε (slope 0.54–0.95), and the chaotic-layer
transition κ⊥ ∝ D_QL (constant ratio ≈ 0.12 across a 16× range of overlap). The
discrete maximum principle is satisfied and solver residuals are ≤ 10⁻⁸, except
the single hardest solve (m = 36 at κ⊥ = 10⁻⁶, ~5×10⁵ unknowns, anisotropy 10⁶),
which reached 10⁻⁶ at the iteration cap; the V_PD/ΔT ordering is well clear of
that uncertainty.

## An error that was caught

The first V_PD run used a magnetic field with **constant perturbation
amplitudes**, omitting Paul's ψ(ψ − ψ̄) envelope (and using a different ι profile
and amplitude normalisation, which also changed the chain counts to 7/11/19 for
m = 12/20/36).

**All validation gates passed on the wrong field.** The gates probe a single
resonance, where the envelope is nearly constant, so they were blind to the
difference. The actual defect was at the boundaries: the angle-averaged χ(ρ) for
m = 12 at κ⊥ = 10⁻⁶ is ≈ 0.5–0.7 at the domain edges with the constant-amplitude
field, versus ≈ 0.00 with Paul's envelope. A nonzero χ at the edge means field
lines pierce the T = 0 and T = 1 boundaries — spurious parallel flux into the
volume, which the envelope exists precisely to forbid.

The error was found by a consistency check against the original C2/C3a code
(detection agreement 100%, regular/chaotic classification agreement 93–96%), and
the run was repeated with the correct field. **The conclusions did not change; the
numbers became cleaner** — on Paul's field the ΔT ordering is strictly monotone at
every κ⊥. The superseded record is retained (`RECORD_C3b_vpd.md`, marked
SUPERSEDED).

The lesson: passing a validation does not mean the model is right. The validation
may be blind to the difference.

## What remains open

**The cantori hypothesis.** Why does converse-KAM fail where V_PD succeeds? One
possibility is that converse-KAM answers only whether a KAM surface exists, and
cantori are not KAM surfaces. The m = 36 chaotic sea is connected but transport
through it is slow: a single orbit needs about 10⁵ toroidal turns to span the full
range [1/8, 7/8], versus a width of 0.19 after 10³ turns. This is a hypothesis,
not a result; testing it requires measuring the cantorus turnstile flux directly.

**Island-centred foliation.** Kallinikos, MacKay & Martínez-del-Río (§5) propose
that a foliation centred on the elliptic field lines of an island chain would let
converse-KAM distinguish islands from chaos. Not attempted here. (In the present
radial-direction test, converse-KAM and WBA disagree precisely on islands —
converse-KAM marks an island interior "non-existent", WBA marks it regular.)

**Finite β.** Every analysis here is on a *given* magnetic field. The problem of
pressure modifying the field requires a finite-β 3D equilibrium — HINT.

## Repository layout

| Path | Contents |
|---|---|
| `tools/` | Field models (KMM, Paul), converse-KAM, WBA, Greene residue |
| `vpd/` | V_PD solver (anisotropic heat diffusion) and Paul-field drivers |
| `consistency/` | Cross-check of reimplemented diagnostics against the originals |
| `records/` | Eight records + open questions (Korean) |
| `results/grid/` | Committed grid data (t_c, dig, χ, T) |
| `figures/` | Poincaré sections, detection maps, dig maps, comparison plots |
| `scripts/` | Analysis and figure-generation drivers |

Grid data is committed, so the three-axis comparison can be reproduced without
rerunning the solvers.

## Reproduce

```bash
git clone https://github.com/kmh21012799-cpu/sv-4
cd sv-4
pip install -r requirements.txt

# Fields and converse-KAM / WBA validation gates
python scripts/poincare_paul.py            # Poincaré sections of the four fields
python scripts/converse_kam_kmm.py         # C1 island-area / residue validation
python scripts/wba_gate2.py                # WBA gate (regular vs chaotic)

# V_PD on Paul's field
python vpd/validate1_paul.py               # single-resonance gate
python vpd/validate2_paul.py               # chaotic-layer gate
python vpd/stage2_paul.py                  # V_PD, ΔT for the four fields
python vpd/stage3_paul.py                  # three-axis table (one field, original diagnostics)
python vpd/correlations_paul.py            # r(converse-KAM, V_PD) vs κ⊥
```

## Related repositories

- [sv1_comp](https://github.com/kmh21012799-cpu/sv1_comp) — full narrative
- [sv-1](https://github.com/kmh21012799-cpu/sv-1) — standard-map diffusion
- [sv-2](https://github.com/kmh21012799-cpu/sv-2) — QUASR vacuum configurations
- [sv-3](https://github.com/kmh21012799-cpu/sv-3) — WBA

## References

- MacKay, R. S. (2018) *Reg. Chaotic Dyn.* **23**, 797
- Kallinikos, N., MacKay, R. S. & Martínez-del-Río, D. (2023) *PPCF* **65**, 095021
- **Paul, E. J., Hudson, S. R. & Helander, P. (2022) *J. Plasma Phys.* 88, 905880107**
- Duignan, N. & Meiss, J. D. (2023) *Physica D* **449**, 133749
- Vlad, M. et al. (2002)

## Provenance

The three diagnostics reproduce their source papers rather than reusing external
code; `tools/SOURCE.txt` records this, and `consistency/` documents the
cross-check of the V_PD-session reimplementations against the C2/C3a originals.
The staged records (C0–C3b) and `future_questions.md` are in Korean.
