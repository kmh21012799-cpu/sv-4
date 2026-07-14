"""
C2-C: converse-KAM on the four Paul critical-overlap fields.
Runs ONE field to t_f (recording t_c per point), saves the map, and reports all
metrics.  Because t_c is recorded, area{t_c<=tf} and timeout fractions for
tf=200/400/800 are cumulative from a single t_f=800 run -- no re-runs.

Area is a flux area in (rho,theta); we report the full domain and Paul's
rho in [0.25,0.75] sub-band.  Primary signal is the t_c distribution (median,
IQR, histogram), which is far less sensitive to boundary aliasing than area.
"""
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
from field_paul import paul_m4, paul_m12, paul_m20, paul_m36        # noqa: E402
from converse_kam import detection_map_uv, nonexistence_area_uv     # noqa: E402
DATA = os.path.join(os.path.dirname(__file__), "..", "data")
FIELDS = {"m4": paul_m4, "m12": paul_m12, "m20": paul_m20, "m36": paul_m36}


def run(key, N=160, t_f=800.0, tag="tc"):
    f = FIELDS[key]()
    rho = np.linspace(0.0, 1.0, N)
    th = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False)   # periodic: no dup endpoint
    t0 = time.time()
    tc = detection_map_uv(f, rho, th, t_f=t_f, symmetry=False)
    dt = time.time() - t0
    np.savez(os.path.join(DATA, f"paul_{key}_{tag}.npz"), tc=tc, rho=rho, th=th, t_f=t_f)
    report(key, tc, rho, th, dt)


def report(key, tc, rho, th, dt=None):
    valid = ~np.isnan(tc)                       # exclude the rho-boundary guard
    nvalid = valid.sum()
    full_dom = 1.0 * 2.0 * np.pi                 # true (rho,theta) domain area
    band = (rho >= 0.25) & (rho <= 0.75)
    line = [f"[{key}] N={len(rho)}"]
    if dt is not None:
        line.append(f"{dt:.0f}s")
    print(" ".join(line))
    for tf in [200.0, 400.0, 800.0]:
        det = valid & np.isfinite(tc) & (tc <= tf)
        undet = valid & ~det
        S_full = nonexistence_area_uv(np.where(det, tc, np.inf), rho, th, tf)
        # sub-band as a clean DETECTED FRACTION (<=100%), robust to grid endpoints
        vb = valid[band, :]; db = det[band, :]
        subfrac = 100.0 * db.sum() / vb.sum()
        tcv = tc[det]
        med = np.median(tcv) if tcv.size else np.nan
        q1, q3 = (np.percentile(tcv, [25, 75]) if tcv.size else (np.nan, np.nan))
        print(f"   t_f={tf:4.0f}: area_full={S_full:.4f} ({100*S_full/full_dom:.1f}%dom)"
              f"  detected[.25,.75]={subfrac:.1f}%"
              f"  undetected={100*undet.sum()/nvalid:.1f}%"
              f"  t_c med={med:.1f} IQR=[{q1:.1f},{q3:.1f}]")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("key")
    ap.add_argument("--N", type=int, default=160)
    ap.add_argument("--tf", type=float, default=800.0)
    ap.add_argument("--tag", default="tc")
    ap.add_argument("--report-only", action="store_true")
    args = ap.parse_args()
    if args.report_only:
        d = np.load(os.path.join(DATA, f"paul_{args.key}_{args.tag}.npz"))
        report(args.key, d["tc"], d["rho"], d["th"])
    else:
        run(args.key, N=args.N, t_f=args.tf, tag=args.tag)
