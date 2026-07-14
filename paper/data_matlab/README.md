# MATLAB data export — corrected-field (Paul envelope) results

MATLAB `.mat` files for paper figure generation, converted from the committed
data with `scipy.io.savemat`. **All data are from the corrected field**
(Paul's Hamiltonian with the `psi(psi - psibar)` envelope). The superseded
constant-amplitude-field data (`results/stage2.json`, `stage3.json`,
`validation1.json`, `validation2.json`) is **not** exported.

Load in MATLAB with `load('stage3_m4.mat')` etc.

## Conventions and caveats (apply throughout)

- **rho** in [0.25, 0.75] — the transport-relevant core, **common to all fields**.
- **theta** in [0, 2*pi/m) — the reduced poloidal cell; **its range differs per
  field** (m = 4: [0, pi/2); m = 36: [0, pi/18)). All modes share the poloidal
  number m, so the solution is periodic with period 2*pi/m. Do not assume a
  common theta axis across fields.
- **chi** is **binary (0 or 1)** — the local V_PD indicator
  `Theta(kpar |grad_par T|^2 - kperp |grad_perp T|^2)`. **V_PD** (the scalar in
  `stage2_paul.mat`) is the fraction of core points with chi = 1.
- **tc** contains **Inf** at undetected points (no converse-KAM detection within
  t_f = 200). Use the **detected** mask (1 = detected) rather than filtering Inf
  by hand. MATLAB understands Inf.
- **T** is a single **zeta = 0 slice**, not a 3-D field.
- Spatial maps of **chi and T exist only at kperp = 1e-6** (the small-kperp,
  most parallel-dominated case). t_c and dig are field-line diagnostics and are
  **kperp-independent**.
- Source repo: `kmh21012799-cpu/sv-4`, branch `main`. Corrected-field data was
  produced in C3b v2 (commit `8a52ca1`); converse-KAM/WBA source code is
  `tools/` (commit `393328e` / `bf5ca06`).

---

## 1. `stage3_m{4,12,20,36}.mat` — three-axis grid maps (the core figure)

One file per field. Same (rho, theta) grid for all three axes.
Source: `results/grid/stage3_m{m}.npz`.

| Variable | Size | Units | Meaning |
|---|---|---|---|
| `tc` | 24 x 24 | zeta (toroidal angle) | converse-KAM first-detection time; **Inf = undetected** |
| `dig` | 24 x 24 | digits (dimensionless) | WBA convergence: high (~10-14) regular, low (~1) chaotic |
| `chi` | 24 x 24 | 0 / 1 | V_PD local indicator (1 = parallel diffusion dominates) at kperp = 1e-6 |
| `detected` | 24 x 24 | 0 / 1 | converse-KAM detection mask (1 = tc finite) |
| `rho` | 24 x 1 | — | radial coordinate, [0.25, 0.75] |
| `theta` | 24 x 1 | rad | poloidal coordinate, [0, 2*pi/m) |

Indexing: `tc(i,j)` at `rho(i)`, `theta(j)`.

## 2. `vpd_m{4,12,20,36}.mat` — high-resolution T / chi slice

One file per field, zeta = 0, kperp = 1e-6. For isotherms and the
parallel-dominant region (cf. Paul Fig. 6). Source:
`results/grid/vpd_m{m}_kperp1e-06.npz`.

| Variable | Size | Units | Meaning |
|---|---|---|---|
| `T_zeta0` | 129 x 32 | normalised (T=0 inner, T=1 outer) | temperature field, zeta = 0 slice |
| `chi_zeta0` | 129 x 32 | 0 / 1 | V_PD indicator on the same slice |
| `rho` | 129 x 1 | — | radial coordinate, [0, 1] (full domain here, not just the core) |
| `theta` | 32 x 1 | rad | poloidal coordinate, [0, 2*pi/m) |

Note: the continuous diffusion ratio itself is not stored — only its sign
(`chi_zeta0`). Reconstructing the ratio would need the 3-D T (a recompute).

## 3. `stage2_paul.mat` — V_PD and Delta T vs kperp

Source: `results/stage2_paul.json`. Rows are fields, columns are kperp.

| Variable | Size | Meaning |
|---|---|---|
| `m` | 4 x 1 | poloidal number [4; 12; 20; 36] |
| `kperp` | 3 x 1 | [1e-4; 1e-5; 1e-6] |
| `V_PD` | 4 x 3 | parallel-diffusion volume fraction; `V_PD(i,j)` = field m(i), kperp(j) |
| `DeltaT` | 4 x 3 | angle-averaged core temperature drop <T(0.75)> - <T(0.25)> |

## 4. `correlations_paul.mat` — pointwise correlations vs kperp

Source: `results/correlations_paul.json`. Rows fields, columns kperp = [1e-4; 1e-5; 1e-6].

| Variable | Size | Meaning |
|---|---|---|
| `m` | 4 x 1 | [4; 12; 20; 36] |
| `kperp` | 3 x 1 | [1e-4; 1e-5; 1e-6] |
| `r_tc_vpd` | 4 x 3 | Spearman r(-tc, chi): converse-KAM vs V_PD (key question) |
| `r_wba_vpd` | 4 x 3 | Spearman r(dig, chi): WBA vs V_PD |
| `r_tc_wba` | 4 x 3 | Spearman r(-tc, dig): converse-KAM vs WBA |

**`r_tc_wba` is kperp-independent** (both axes are field-line diagnostics); it is
tiled across the 3 columns, so the columns are identical. `r_tc_vpd` does not grow
toward small kperp for any field — this is the refutation of Paul's expectation.

## 5. `validation1_paul.mat` — single-island gate (Paul Fig. 1)

Source: `results/validation1_paul.json`. V_PD vs perturbation amplitude eps for a
single (m=2, n=1) island, three kperp.

| Variable | Size | Meaning |
|---|---|---|
| `eps` | 16 x 3 | amplitude scan. **NOT a shared 16 x 1 vector**: the grid differs per kperp (each is log-spaced around sqrt(kperp)), so `eps(:,j)` goes with `V_PD(:,j)` |
| `kperp` | 3 x 1 | **[1e-2; 1e-4; 1e-6]** (this file's kperp set, not 1e-5) |
| `V_PD` | 16 x 3 | V_PD at (eps(i,j), kperp(j)) |
| `eps_crit` | 3 x 1 | sqrt(kperp)/2, the predicted transition |

## 6. `validation2_paul.mat` — chaotic-layer gate (Paul Fig. 3)

Source: `results/validation2_paul.json`. V_PD vs kperp for a chaotic layer (m=12,
n=2..10) at three Chirikov overlap parameters.

| Variable | Size | Meaning |
|---|---|---|
| `kperp` | 9 x 1 | scan, 1e-2 ... 1e-6 |
| `S` | 3 x 1 | Chirikov overlap [2; 2.828; 4] |
| `V_PD` | 9 x 3 | `V_PD(i,j)` = kperp(i), S(j) |
| `D_QL` | 3 x 1 | quasilinear field-line diffusion coefficient per S |
| `kperp_trans` | 3 x 1 | kperp where V_PD crosses 0.5, per S |

## 7. `wba_tconv.mat` — WBA dig T-convergence (curve crossing)

Source: **`records/RECORD_C3a_wba_paul.md`** table (commit `bf5ca06`), read
verbatim — the per-point data was not saved, only these medians.

| Variable | Size | Meaning |
|---|---|---|
| `m` | 4 x 1 | [4; 12; 20; 36] |
| `T` | 4 x 1 | integration time [500; 1000; 2000; 5000] (toroidal periods) |
| `dig_median` | 4 x 4 | median core dig; `dig_median(i,j)` = field m(i), time T(j) |

m=36 falls (1.42 -> 0.99) while m=4 rises (1.05 -> 1.46) as T grows: the ordering
at short T is a finite-time artifact and the curves cross.

## 8. `poincare_m{4,12,20,36}.mat` — Poincare sections (recomputed)

Section zeta = 0, on the corrected envelope field. **Recomputed** with
`paper/make_poincare.py` (the orbit points were never saved). Arrays also live in
`paper/data/poincare_m{m}.npz` (committed for reproducibility). m=20 was newly
computed (it had no prior figure).

| Variable | Size | Meaning |
|---|---|---|
| `rho` | Npts x 1 | psi = rho in [0, 1] (full domain) |
| `theta` | Npts x 1 | poloidal angle in [0, 2*pi) (physical, not reduced) |
| `kind` | Npts x 1 | 0 = regular seed (theta=0 line), 1 = chaos-filling orbit |
| `res_psi` | (#chains) x 1 | resonance locations n/m (for reference lines) |

Npts differs per field. Note theta here spans the full [0, 2*pi) (Poincare of the
full map), unlike the reduced-cell maps in files 1-2.
