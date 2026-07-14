"""
Consistency check: ORIGINAL (C2/C3a) converse-KAM & WBA  vs  C3b reimplementation,
on the SAME C3b field, SAME (rho,theta) grid points.

Isolates the DIAGNOSTIC-ALGORITHM difference (field is held identical).
Outputs Spearman, Pearson, detection/chaos agreement, and absolute-scale ratio.
"""
import sys, os, json, time
import numpy as np
from multiprocessing import Pool
from scipy.stats import spearmanr, pearsonr

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "orig"))
import converse_kam as orig_ck          # noqa: E402
import wba as orig_wba                  # noqa: E402

from vpd.field import paul_field        # noqa: E402
from vpd.diagnostics import run_diagnostics  # noqa: E402
from consistency.adapter import C3bAdapter    # noqa: E402


def _ck_worker(args):
    field, rho0, th0, t_f = args
    return orig_ck.detect_tc(field, rho0, th0, t_f=t_f, symmetry=False)


def _wba_worker(args):
    field, rho0, th0, n_periods = args
    return orig_wba.wba_point(field, rho0, th0, n_periods=n_periods)["dig_psi"]


def run(m=12, n_rho=28, n_theta=24, t_f=200.0, wba_periods=400, nproc=4):
    f = paul_field(m, chirikov=1.0)
    adap = C3bAdapter(f)
    rho = np.linspace(0.25, 0.75, n_rho)
    theta = np.linspace(0, 2 * np.pi / m, n_theta, endpoint=False)
    RHO, TH = np.meshgrid(rho, theta, indexing="ij")
    pts = list(zip(RHO.ravel(), TH.ravel()))

    # ---- ORIGINAL converse-KAM t_c (in zeta/phi units) ----
    t0 = time.time()
    with Pool(nproc) as p:
        tc_orig = np.array(p.map(_ck_worker, [(adap, r, t, t_f) for r, t in pts]))
    print(f"orig converse-KAM: {time.time()-t0:.0f}s")

    # ---- ORIGINAL WBA dig ----
    t0 = time.time()
    with Pool(nproc) as p:
        dig_orig = np.array(p.map(_wba_worker, [(adap, r, t, wba_periods) for r, t in pts]))
    print(f"orig WBA: {time.time()-t0:.0f}s")

    # ---- C3b reimplementation on the same grid ----
    t0 = time.time()
    d = run_diagnostics(f, rho, theta, n_periods=wba_periods, substeps=28,
                        ck_periods=int(t_f / (2 * np.pi)) + 2)
    print(f"C3b diagnostics: {time.time()-t0:.0f}s")
    tc_mine = d["t_c"].ravel() * (2 * np.pi)   # periods -> zeta, match orig units
    dig_mine = d["dig"].ravel()
    det_mine = d["detected"].ravel()

    det_orig = np.isfinite(tc_orig) & (tc_orig <= t_f)

    out = {"m": m, "n_points": len(pts), "t_f": t_f, "wba_periods": wba_periods}

    # converse-KAM comparison (on points detected by BOTH, for rank/scale)
    both = det_orig & (det_mine > 0.5) & np.isfinite(tc_mine)
    out["ck_detect_frac_orig"] = float(det_orig.mean())
    out["ck_detect_frac_mine"] = float(det_mine.mean())
    out["ck_detect_agreement"] = float(np.mean(det_orig == (det_mine > 0.5)))
    if both.sum() > 5:
        out["ck_spearman"] = float(spearmanr(tc_orig[both], tc_mine[both])[0])
        out["ck_pearson"] = float(pearsonr(tc_orig[both], tc_mine[both])[0])
        out["ck_scale_ratio_median"] = float(np.median(tc_mine[both] / tc_orig[both]))
        out["ck_tc_median_orig"] = float(np.median(tc_orig[det_orig]))
        out["ck_tc_median_mine"] = float(np.median(tc_mine[det_mine > 0.5]))
    # WBA comparison (all points)
    finite = np.isfinite(dig_orig) & np.isfinite(dig_mine)
    out["wba_spearman"] = float(spearmanr(dig_orig[finite], dig_mine[finite])[0])
    out["wba_pearson"] = float(pearsonr(dig_orig[finite], dig_mine[finite])[0])
    out["wba_chaos_frac_orig"] = float(np.mean(dig_orig[finite] < 5))
    out["wba_chaos_frac_mine"] = float(np.mean(dig_mine[finite] < 5))
    out["wba_dig_median_orig"] = float(np.median(dig_orig[finite]))
    out["wba_dig_median_mine"] = float(np.median(dig_mine[finite]))
    out["wba_regular_frac_orig"] = float(np.mean(dig_orig[finite] > 8))
    out["wba_regular_frac_mine"] = float(np.mean(dig_mine[finite] > 8))

    # save raw arrays for a scatter figure
    np.savez("results/consistency_m%d.npz" % m, tc_orig=tc_orig, tc_mine=tc_mine,
             dig_orig=dig_orig, dig_mine=dig_mine, det_orig=det_orig,
             det_mine=det_mine, rho=rho, theta=theta)
    return out


if __name__ == "__main__":
    m = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    res = run(m=m)
    with open(f"results/consistency_m{m}.json", "w") as fh:
        json.dump(res, fh, indent=2)
    print(json.dumps(res, indent=2))
