"""How does the ideal N prediction improve with a longer record? (contagious world)"""
import json, numpy as np
from scipy.optimize import least_squares
import reference_n as R
from observe import lonlat, visible, CLEAR, NOISE
tr = np.load("../private/truth_run_contagious.npz", allow_pickle=True); names = list(tr["names"]); iN = names.index("N")
rng = np.random.default_rng(2026); clear_night = rng.random(43831) < CLEAR
vis = visible(tr["X"], 2, 0, iN, True) & clear_night
tlo, tla = lonlat(tr["X"][:, iN] - tr["X"][:, 2])
truth = {d: (lo, la) for d, b, lo, la in json.load(open("../private/truth_predictions_contagious.json"))["predictions_3.csv"]}
pdays = np.array(sorted(truth)); Tlo = np.array([truth[d][0] for d in pdays]); Tla = np.array([truth[d][1] for d in pdays])
WINDOWS = {"last 60 nights": [(36466, 36525)], "current appearance (yr 99.2-100)": [(36233, 36525)],
           "current + previous (yr 71-72)": [(26050, 26343), (36233, 36525)],
           "all four appearances": [(5666, 5960), (15836, 16130), (26050, 26343), (36233, 36525)]}
for wname, wins in WINDOWS.items():
    days = np.concatenate([np.arange(a, b + 1) for a, b in wins]); days = days[vis[days]]
    nrng = np.random.default_rng(5)
    olat = tla[days] + nrng.normal(0, np.degrees(NOISE), len(days))
    olon = tlo[days] + nrng.normal(0, np.degrees(NOISE), len(days)) / np.cos(np.radians(tla[days]))
    E = int(days[0])
    R.EPOCH = E; R.B_X, R.B_V, R.B_Q = tr["X"][E, :8], tr["V"][E, :8], tr["Q"][E, :8]
    x_true = np.r_[tr["X"][E, iN], tr["V"][E, iN]]; f_true = tr["Q"][E, iN] / tr["m"][iN]
    rel = days - E
    def fitlaw(free_f):
        def resid(p):
            r = R.propagate(p[:6][None], np.array([p[6] if free_f else 0.0]), 0.0, 36525 - E, sample_days=1)
            lo, la = R.directions(r, 1); lo, la = lo[rel, 0], la[rel, 0]
            return np.r_[((lo - olon + 180) % 360 - 180) * np.cos(np.radians(olat)), la - olat] * 60
        p0 = np.r_[x_true, np.log(f_true)] if free_f else x_true
        s = least_squares(resid, p0, x_scale=np.r_[np.full(3, 1e-3), np.full(3, 1e-2), [1e-2]][:len(p0)], diff_step=1e-7,
                          xtol=1e-12, ftol=1e-12, gtol=1e-12)
        cov = np.linalg.pinv(s.jac.T @ s.jac) * (s.fun ** 2).sum() / (len(s.fun) - len(p0))
        return s.x, cov, np.sqrt((s.fun ** 2).mean())
    line = [f"{wname}: {len(days)} sightings"]
    for ref, lam, free_f in [("true law", 50.0, True), ("no exchange", 0.0, True), ("plain Newton", 0.0, False)]:
        p, cov, rms = fitlaw(free_f)
        r = R.propagate(p[:6][None], np.array([p[6] if free_f else 0.0]), lam, 43830 - E, sample_days=1)
        lo, la = R.directions(r, 1); idx = pdays - E
        err = R.angsep(lo[idx, 0], la[idx, 0], Tlo, Tla)
        m = {k: np.median(err[np.abs(pdays - d) <= 15]) for k, d in R.MARKS.items()}
        pf = f"pull {np.exp(p[6]):.3f}+/-{np.exp(p[6])*np.sqrt(cov[6,6]):.3f}" if free_f else "pull 1 (fixed)"
        line.append(f"   {ref:12s} {pf:22s} rms {rms:.2f}'  errors: 1 yr {m['1 yr']:.3f}, 5 yr {m['5 yr']:.3f}, 20 yr {m['20 yr']:.3f} deg")
    print("\n".join(line))
