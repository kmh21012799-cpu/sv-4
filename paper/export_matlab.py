"""
Export the CORRECTED (Paul-envelope field) committed data to MATLAB .mat files
for figure generation. Conversion only -- no new physics, no plotting.

Sources (all committed, corrected field):
  results/grid/stage3_m*.npz, results/grid/vpd_m*_kperp1e-06.npz,
  results/stage2_paul.json, results/correlations_paul.json,
  results/validation1_paul.json, results/validation2_paul.json,
  records/RECORD_C3a_wba_paul.md (T-convergence table, read verbatim).

SUPERSEDED (old constant-amplitude field) data is NOT exported.
"""
import json
import numpy as np
from scipy.io import savemat

OUT = "paper/data_matlab"
MS = [4, 12, 20, 36]
KPERP3 = [1e-4, 1e-5, 1e-6]          # stage2 / correlations column order


# ---- 1. three-axis grid maps: stage3_m{m}.mat ---------------------------
for m in MS:
    d = np.load(f"results/grid/stage3_m{m}.npz")
    savemat(f"{OUT}/stage3_m{m}.mat", {
        "tc": d["tc"].astype(float),          # Inf at undetected points (kept)
        "dig": d["dig"].astype(float),
        "chi": d["chi"].astype(float),        # binary 0/1
        "detected": d["detected"].astype(float),
        "rho": d["rho"].reshape(-1, 1),
        "theta": d["theta"].reshape(-1, 1),
    }, do_compression=True)
print("wrote stage3_m{4,12,20,36}.mat")

# ---- 2. high-res T / chi slice: vpd_m{m}.mat ----------------------------
for m in MS:
    d = np.load(f"results/grid/vpd_m{m}_kperp1e-06.npz")
    savemat(f"{OUT}/vpd_m{m}.mat", {
        "T_zeta0": d["T_zeta0"].astype(float),
        "chi_zeta0": d["chi_zeta0"].astype(float),
        "rho": d["rho"].reshape(-1, 1),
        "theta": d["theta"].reshape(-1, 1),
    }, do_compression=True)
print("wrote vpd_m{4,12,20,36}.mat")

# ---- 3. stage2_paul: V_PD, DeltaT (4 x 3) -------------------------------
s2 = json.load(open("results/stage2_paul.json"))
VPD = np.array([[next(r for r in s2[str(m)]["rows"] if abs(r["kperp"] - kp) < 1e-12)["V_PD"]
                 for kp in KPERP3] for m in MS])
DT = np.array([[next(r for r in s2[str(m)]["rows"] if abs(r["kperp"] - kp) < 1e-12)["DeltaT"]
                for kp in KPERP3] for m in MS])
savemat(f"{OUT}/stage2_paul.mat", {
    "m": np.array(MS, float).reshape(-1, 1),
    "kperp": np.array(KPERP3, float).reshape(-1, 1),
    "V_PD": VPD, "DeltaT": DT,
})
print("wrote stage2_paul.mat  V_PD/DeltaT shape", VPD.shape)

# ---- 4. correlations_paul: r_tc_vpd, r_wba_vpd (4x3), r_tc_wba (4x1->4x3) --
cr = json.load(open("results/correlations_paul.json"))
kk = ["1e-04", "1e-05", "1e-06"]
r_tc_vpd = np.array([[cr[str(m)]["per_kperp"][k]["r_tc_vpd"] for k in kk] for m in MS])
r_wba_vpd = np.array([[cr[str(m)]["per_kperp"][k]["r_wba_vpd"] for k in kk] for m in MS])
r_tc_wba_col = np.array([cr[str(m)]["r_tc_wba"] for m in MS]).reshape(-1, 1)  # kperp-independent
savemat(f"{OUT}/correlations_paul.mat", {
    "m": np.array(MS, float).reshape(-1, 1),
    "kperp": np.array(KPERP3, float).reshape(-1, 1),
    "r_tc_vpd": r_tc_vpd,
    "r_wba_vpd": r_wba_vpd,
    "r_tc_wba": np.tile(r_tc_wba_col, (1, 3)),   # tiled to 4x3 (columns identical)
})
print("wrote correlations_paul.mat")

# ---- 5. validation1_paul: eps (16x3, differs per kperp), V_PD (16x3) ----
v1 = json.load(open("results/validation1_paul.json"))
v1keys = ["1e-02", "1e-04", "1e-06"]             # this file's kperp keys
eps16 = np.array([[r["eps"] for r in v1[k]["rows"]] for k in v1keys]).T   # 16 x 3
vpd16 = np.array([[r["V_PD"] for r in v1[k]["rows"]] for k in v1keys]).T
savemat(f"{OUT}/validation1_paul.mat", {
    "eps": eps16,                                 # 16 x 3 (per-kperp grids)
    "kperp": np.array([1e-2, 1e-4, 1e-6]).reshape(-1, 1),
    "V_PD": vpd16,
    "eps_crit": np.array([v1[k]["eps_crit"] for k in v1keys]).reshape(-1, 1),
})
print("wrote validation1_paul.mat  eps/V_PD shape", eps16.shape)

# ---- 6. validation2_paul: kperp(9), S(3), V_PD(9x3), D_QL(3), trans(3) ---
v2 = json.load(open("results/validation2_paul.json"))
Skeys = ["2.000", "2.828", "4.000"]
kperp9 = np.array([r["kperp"] for r in v2[Skeys[0]]["rows"]])
vpd93 = np.array([[r["V_PD"] for r in v2[s]["rows"]] for s in Skeys]).T     # 9 x 3
savemat(f"{OUT}/validation2_paul.mat", {
    "kperp": kperp9.reshape(-1, 1),
    "S": np.array([v2[s]["chirikov"] for s in Skeys]).reshape(-1, 1),
    "V_PD": vpd93,
    "D_QL": np.array([v2[s]["D_QL"] for s in Skeys]).reshape(-1, 1),
    "kperp_trans": np.array([v2[s]["kperp_trans"] for s in Skeys]).reshape(-1, 1),
})
print("wrote validation2_paul.mat  V_PD shape", vpd93.shape)

# ---- 7. WBA T-convergence (RECORD_C3a table, read verbatim) -------------
# | field | T=500 | T=1000 | T=2000 | T=5000 |
# | m=4   | 1.05  | 1.19   | 1.21   | 1.46   |
# | m=12  | 0.93  | 0.93   | 0.83   | 0.97   |
# | m=20  | 1.21  | 1.12   | 1.01   | 0.93   |
# | m=36  | 1.42  | 1.24   | 1.14   | 0.99   |
dig_median = np.array([
    [1.05, 1.19, 1.21, 1.46],
    [0.93, 0.93, 0.83, 0.97],
    [1.21, 1.12, 1.01, 0.93],
    [1.42, 1.24, 1.14, 0.99],
])
savemat(f"{OUT}/wba_tconv.mat", {
    "m": np.array(MS, float).reshape(-1, 1),
    "T": np.array([500, 1000, 2000, 5000], float).reshape(-1, 1),
    "dig_median": dig_median,                     # 4 (fields) x 4 (T)
})
print("wrote wba_tconv.mat")
print("ALL CONVERSIONS DONE")
