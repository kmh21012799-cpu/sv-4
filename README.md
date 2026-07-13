# converse-kam-3d

Reproduction/verification project working toward a quantitative comparison of
**Paul–Hudson–Helander's effective parallel-diffusion volume `V_PD`** with the
**converse-KAM** non-existence method, in the *same* magnetic fields.

Paul, Hudson & Helander (2022, arXiv:2108.06328) wrote that their metric
"appears to be related to the converse KAM approach … We reserve such a
comparison for future publication." This repo builds toward that comparison in
staged, verifiable steps. **It makes no new-discovery claims.**

## Stages

- **C0 — magnetic fields** (`records/RECORD_C0_fields.md`): implement and
  Poincaré-validate the KMM (2023) and Paul (2022) model fields. **Done / PASS.**
- **C1 — converse-KAM 3D** (`records/RECORD_C1_converse_kam.md`): implement the
  MacKay/KMM converse-KAM test and validate it on the integrable KMM island
  (known answer).
- C2 / C3 — apply to Paul's critical-chaos fields; implement `V_PD` and compare.
  **Out of scope for now.**

## Layout

```
tools/     field_kmm.py, field_paul.py, kmm_island.py, residue.py, converse_kam.py
scripts/   Poincaré / residue / converse-KAM figure generators
figures/   generated PNGs
records/   RECORD_C0_fields.md, RECORD_C1_converse_kam.md, future_questions.md
```

## Setup

```
pip install -r requirements.txt
```

See each `RECORD_*.md` for exact reproduction commands and the LIMITATIONS
(stated first in every record).
