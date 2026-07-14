# RECORD A2 — field-model discrepancy (critical)

**Character: verification finding. Code inspection only, no new computation.
Neutral, honest. Verdict: CASE 1.**

A final check of the three-axis comparison asked whether C2/C3a and C3b used the
*same* magnetic field. Paul's specification is a single one (paper eq. 4.1):

    chi = iota' psi^2/2  +  sum_mn eps_mn * psi(psi - psibar) * cos(m theta - n zeta)

with the `psi(psi - psibar)` **envelope**, chosen so the perturbation vanishes
at both boundaries (rho=0 and rho=rhobar) — "no parallel flux into the volume".

## Verdict: CASE 1 — C3b used a genuinely different field (envelope missing).

CASE 2 (mere notation) is ruled out by the code. From
`consistency/orig/field_paul.py` (C2/C3a) vs `vpd/field.py` (C3b):

| feature | ORIGINAL C2/C3a | C3b |
|---|---|---|
| `psi(psi-psibar)` envelope | **present**: `dpsi = Σε m sin · psi(psi-psibar)` | **absent**: `drho = Σε m sin` |
| `dtheta/dzeta` perturbation | present: `(2psi-psibar)Σε cos` | absent: `dtheta = iota(rho)` only |
| iota profile | `iota = psi` (0 at psi=0) | `iota = 0.5 + rho` |
| critical amplitude | `1/(16 n(m-n))` (n-dependent) | `1/(16 m²)` (uniform; ×4 smaller at mid-n) |
| # island chains (m=12/20/36) | 9 / 15 / 27 | **7 / 11 / 19** |
| resonance radial span (m=36) | [0.14, 0.86] | [0.25, 0.75] only |

(m=4 alone happens to share resonance *locations* [0.25,0.5,0.75], but still
differs in amplitude ×4, envelope, and the dtheta term.)

### Why this is CASE 1, not CASE 3 (harmless)
1. **Boundary parallel flux.** C3b's perturbation is nonzero at rho=0,1, so
   `drho/dzeta ≠ 0` there and field lines pierce the V_PD Dirichlet boundaries
   (T=0, T=1). Paul's envelope exists precisely to forbid this. This is a defect
   in the V_PD computation itself, not a cosmetic difference.
2. **Number of island chains differs** (m=12/20/36: 7/11/19 vs 9/15/27). The
   chain count is the very variable Paul's m-scan probes; it is not the same
   experiment.
3. **Envelope varies ~33% across the core** (`psi(psi-1)`: −0.1875 at rho=0.25,0.75
   to −0.25 at rho=0.5 → ~15% island-width variation), on top of the ×4 amplitude
   and iota-offset differences.

## What still holds
- **The V_PD solver is correct.** Stage-1 validation (eps_crit=√kperp/2, √eps
  scaling, chaotic-layer transition) is *local* single-resonance physics,
  insensitive to the global envelope — so it validated the *solver*, and passing
  it did **not** certify the field spectrum (as the check anticipated).
- **RECORD_A (diagnostic consistency) is unaffected**: it explicitly compared
  original vs C3b diagnostics *on the C3b field*, and that comparison stands.
- **The three-axis table is internally self-consistent**: all three indicators
  in C3b were computed on the same C3b field, so it is a valid comparison *on
  that field* — not "meaningless", but on the wrong field.

## What does NOT hold (must be corrected)
- The claim that Stage 2 **reproduces Paul's four fields / Fig 5–6**: the field
  is not Paul's.
- The Stage-2/3 *quantitative* results (V_PD largest at m=4; ΔT largest at m=36;
  the correlations; key question #8) are established **only on the C3b field**;
  whether they hold on Paul's actual field is **unknown**.
- "Completing Paul's deferred comparison **on Paul's fields**": it was done on a
  different field.
- The original C2/C3a results (Paul's field) **cannot** be dropped into the C3b
  table — the table used C3b-field recomputations, not those originals.

## Required correction (NOT executed here — reported and stopped per instruction)
Recompute V_PD (re-run Stage-1 gate, Stage 2, Stage 3) on Paul's actual field
(`consistency/orig/field_paul.py`, envelope included). On that field the
comparison can additionally use C2/C3a's **original** converse-KAM/WBA results
directly, making it a true "three campaigns, one stage". This is a substantial
recomputation and is deferred to an explicit go-ahead.

**This is logged rather than hidden. Redoing is not the shameful thing; hiding
would be.**
