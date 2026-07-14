# sv-4 — C3b: V_PD and the three-axis comparison

Completes the converse-KAM vs V_PD comparison that Paul, Hudson & Helander
(2022) deferred to "future publication", putting three non-integrability
diagnostics on one stage — same fields, same grid, same core:

- **converse-KAM** (topology) — `vpd/diagnostics.py`
- **WBA** (dynamics) — `vpd/diagnostics.py`
- **V_PD**, effective volume of parallel diffusion (transport) — `vpd/solver3d.py`

**Read `RECORD_C3b_vpd.md` first (LIMITATIONS up front).**

## Layout
- `vpd/field.py` — reduced slab field-line model + the four Paul fields
- `vpd/solver2d.py`, `vpd/solver3d.py` — anisotropic heat-diffusion solvers (Q1
  FEM, SPD, AMG-preconditioned CG)
- `vpd/diagnostics.py` — converse-KAM (t_c) and WBA (dig) field-line diagnostics
- `vpd/validate1.py`, `vpd/validate2.py` — Stage-1 gate (Paul Fig 1, Fig 3)
- `vpd/stage2.py` — four fields, V_PD and ΔT (Paul Fig 5, 6)
- `vpd/stage3.py` — the three-axis comparison and correlations
- `vpd/plots.py` — figures; `results/*.json` — data; `figures/*.png` — figures

## Reproduce
```
pip install numpy scipy matplotlib pyamg
python3 -m vpd.validate1     # Val 1 (single island)
python3 -m vpd.validate2     # Val 2 (chaotic layer)
python3 -m vpd.stage2        # four fields
python3 -m vpd.stage3        # three-axis comparison
python3 -m vpd.plots all     # all figures
```

## Headline (neutral)
Transport (V_PD, ΔT) separates Paul's four critical-overlap fields; topology
(converse-KAM) and dynamics (WBA) do not. The converse-KAM↔V_PD correlation does
**not** grow as κ⊥→0, contrary to Paul's stated expectation. Reproduction +
comparison; no new-discovery claim.
