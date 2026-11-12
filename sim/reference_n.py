"""Reference answers for stage 3: how well can N be predicted from its sightings (year 99.2-100), under three laws?

  R_true   : the world's own law (contagious: lam = 50; control: lam = 0), pull factor fitted
  R_static : N keeps its own fitted pull factor, no charge exchange (lam = 0)  [= R_true in the control world]
  R_newton : plain Newton, pull factor fixed at 1, no exchange
  R_exchange (control world only): the charge rule applied where it does not hold

Idealised on purpose (an upper bound on what the data allow): the fitter knows the true positions of S, the
planets, Tellus and the moon on the first night N is seen, and starts its search from N's true state. Each fit
solves for N's position and velocity at the first night of the window (+ log pull factor when free) by least
squares on the angular residuals. Uncertainty: 50 orbits drawn from the fit's covariance, propagated as well.
Score: great-circle error against N's true direction, median over the days within 15 days of 1, 5, 10 and 20 years
after the last sighting (day 36525).
Writes ../private/reference_N.json.
"""
import csv, json, os, numpy as np
from scipy.optimize import least_squares
from world import run
from observe import lonlat

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
TMP = os.path.join(os.path.dirname(os.path.abspath(__file__)), "runs", "tmp_refN")
EPOCH, LAST, END = 36233, 36525, 43830     # EPOCH = first night N is seen
MARKS = {"1 yr": LAST + 365, "5 yr": LAST + 5 * 365, "10 yr": LAST + 10 * 365, "20 yr": END - 15}
key = json.load(open(os.path.join(ROOT, "private", "world_key.json")))["world_key"]
base = np.load(os.path.join(ROOT, "private", "truth_run_contagious.npz"), allow_pickle=True)
bn = list(base["names"]); NB = 8                       # S, P1, T, M, P2, P3, P4, P5
assert bn[:NB] == ["S", "P1", "T", "M", "P2", "P3", "P4", "P5"]
iT = 2
B_X, B_V, B_Q, B_m = base["X"][EPOCH, :NB], base["V"][EPOCH, :NB], base["Q"][EPOCH, :NB], base["m"][:NB]


def propagate(states, logf, lam, days, sample_days=1):
    """Propagate test bodies (each a 6-vector at EPOCH, with pull factor exp(logf)) for `days` days."""
    k = len(states)
    names = bn[:NB] + [f"x{i}" for i in range(k)]
    m = np.r_[B_m, np.full(k, 1e-16)]; q = np.r_[B_Q, np.exp(logf) * 1e-16]
    X = np.vstack([B_X, states[:, :3]]); V = np.vstack([B_V, states[:, 3:]])
    tp = np.r_[np.zeros(NB, int), np.ones(k, int)]
    r = run(names, m, q, X, V, lam=lam, d=0.1, years=days / 365.25, sample_days=sample_days, rtol=1e-12, tp=tp,
            tag="refn", outdir=TMP)
    return r


def directions(r, k):
    d = r["X"][:, NB:NB + k] - r["X"][:, iT][:, None]
    return lonlat(d)                                   # (samples, k) each


def angsep(lo1, la1, lo2, la2):
    lo1, la1, lo2, la2 = map(np.radians, (lo1, la1, lo2, la2))
    c = np.sin(la1) * np.sin(la2) + np.cos(la1) * np.cos(la2) * np.cos(lo1 - lo2)
    return np.degrees(np.arccos(np.clip(c, -1, 1)))


def fit(obs_day, obs_lon, obs_lat, x_true, f_true, lam, free_f):
    rel = obs_day - EPOCH
    def resid(p):
        st = p[:6][None]; lf = p[6] if free_f else 0.0
        r = propagate(st, np.array([lf]), lam, LAST - EPOCH)
        lo, la = directions(r, 1); lo, la = lo[rel, 0], la[rel, 0]
        dlon = ((lo - obs_lon + 180) % 360 - 180) * np.cos(np.radians(obs_lat))
        return np.r_[dlon, la - obs_lat] * 60.0          # arcminutes
    p0 = np.r_[x_true, np.log(f_true)] if free_f else x_true.copy()
    sol = least_squares(resid, p0, x_scale=np.r_[np.full(3, 1e-3), np.full(3, 1e-2), [1e-2]][:len(p0)],
                        diff_step=1e-7, xtol=1e-12, ftol=1e-12, gtol=1e-12)
    J = sol.jac; dof = len(sol.fun) - len(p0)
    s2 = (sol.fun ** 2).sum() / dof
    cov = np.linalg.pinv(J.T @ J) * s2
    return sol.x, cov, np.sqrt(s2)


def main():
    os.makedirs(TMP, exist_ok=True)
    out = {}
    for label, world in key.items():
        tr = np.load(os.path.join(ROOT, "private", f"truth_run_{world}.npz"), allow_pickle=True)
        names = list(tr["names"]); iN = names.index("N")
        x_true = np.r_[tr["X"][EPOCH, iN], tr["V"][EPOCH, iN]]; f_true = tr["Q"][EPOCH, iN] / tr["m"][iN]
        rows = [r for r in csv.DictReader(open(os.path.join(ROOT, "handover", label, "stage3", "observations_3.csv")))]
        od = np.array([int(r["day"]) for r in rows]); olo = np.array([float(r["lon"]) for r in rows])
        ola = np.array([float(r["lat"]) for r in rows])
        truth = {d: (lo, la) for d, b, lo, la in
                 json.load(open(os.path.join(ROOT, "private", f"truth_predictions_{world}.json")))["predictions_3.csv"]}
        pdays = np.array(sorted(truth)); tlo = np.array([truth[d][0] for d in pdays]); tla = np.array([truth[d][1] for d in pdays])
        lam_true = 50.0 if world == "contagious" else 0.0
        res = {}
        refs = [("R_true", lam_true, True), ("R_static", 0.0, True), ("R_newton", 0.0, False)]
        if world == "newton":
            refs.append(("R_exchange", 50.0, True))      # wrongly applying the charge rule in the control world
        for ref, lam, free_f in refs:
            p, cov, rms = fit(od, olo, ola, x_true, f_true, lam, free_f)
            rng = np.random.default_rng(1)
            draws = rng.multivariate_normal(p, cov, 50)
            allp = np.vstack([p, draws])
            st = allp[:, :6]; lf = allp[:, 6] if free_f else np.zeros(len(allp))
            r = propagate(st, lf, lam, END - EPOCH)
            lo, la = directions(r, len(allp))
            idx = pdays - EPOCH
            err = angsep(lo[idx], la[idx], tlo[:, None], tla[:, None])     # (pred days, 1 + 50)
            marks = {}
            for name, md in MARKS.items():
                w = np.abs(pdays - md) <= 15
                e_best = float(np.median(err[w, 0])); e_cloud = np.median(err[w, 1:], axis=0)
                marks[name] = dict(best_fit=round(e_best, 3), cloud_median=round(float(np.median(e_cloud)), 3),
                                   cloud_90pct=round(float(np.percentile(e_cloud, 90)), 3))
            res[ref] = dict(fit_rms_arcmin=round(rms, 3), fitted_pull_factor=round(float(np.exp(p[6])), 4) if free_f else 1.0,
                            pull_factor_sigma=round(float(np.exp(p[6]) * np.sqrt(cov[6, 6])), 4) if free_f else 0.0,
                            errors_deg=marks)
            print(f"{label} ({world}) {ref}: rms {rms:.2f}', pull factor {res[ref]['fitted_pull_factor']} "
                  f"+/- {res[ref]['pull_factor_sigma']} | " +
                  " | ".join(f"{k}: {v['best_fit']:.3f} deg (cloud {v['cloud_median']:.3f}, 90% {v['cloud_90pct']:.3f})"
                             for k, v in marks.items()))
        out[label] = dict(world=world, n_sightings=len(od), true_pull_factor=round(float(f_true), 4), refs=res)
    json.dump(out, open(os.path.join(ROOT, "private", "reference_N.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
