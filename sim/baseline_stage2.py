"""Stage 2 baselines: predict each comet C1-C12 in years 101-120 from its own 100 years of sightings.

  B2_true   : the world's own law (contagious: charge rule, lam = 50; control: lam = 0), pull factor fitted
  B2_static : each comet keeps its own fitted pull factor, no exchange
  B2_newton : pull factor fixed at 1, no exchange
Idealised as for N: true planet, moon and star positions known; search starts at the comet's true state at its
first sighting. Each comet is fitted on its own as a test body (comets do not see each other here, so C11's
meeting with C7 in year 90.7 is not modelled; noted). Parameters: position and velocity at the first sighting
(+ log pull factor when free). Levenberg-Marquardt; the finite-difference Jacobian is computed in one integration
by carrying 7 perturbed copies as test bodies. Uncertainty: 20 draws from the fit's covariance.
Score: great-circle error on the stage 2 prediction days (every 30 days, years 101-120): median per comet.
Writes ../private/baseline_stage2.json.
"""
import json, os, numpy as np
from concurrent.futures import ProcessPoolExecutor
from baseline_lib import ROOT, LAST, END, KEY, truth_run, load_obs, load_truth_pred, integrate, resid_angles, angsep
from observe import lonlat

NB = 8   # S, P1, T, M, P2, P3, P4, P5 come first in the truth run
COMETS = [f"C{k}" for k in range(1, 13)]


def propagate(tr, epoch, P, lam, days):
    """P: (k, 7) rows of [x y z vx vy vz logf]; returns directions from Tellus (samples, k) for days epoch..epoch+days."""
    k = len(P); names = tr["names"][:NB] + [f"x{i}" for i in range(k)]
    m = np.r_[tr["m"][:NB], np.full(k, 1e-16)]; q = np.r_[tr["Q"][epoch, :NB], np.exp(P[:, 6]) * 1e-16]
    X = np.vstack([tr["X"][epoch, :NB], P[:, :3]]); V = np.vstack([tr["V"][epoch, :NB], P[:, 3:6]])
    tp = np.r_[np.zeros(NB, int), np.ones(k, int)]
    r = integrate(names, m, q, X, V, lam, days, tp=tp)
    return lonlat(r["X"][:, NB:] - r["X"][:, 2][:, None])


def fit_comet(args):
    world, label, c, law = args
    tr = truth_run(world); obs = load_obs(label, 2)[c]
    d, olo, ola = obs; epoch = int(d[0]); rel = d - epoch
    ic = tr["names"].index(c)
    lam_true = 50.0 if world == "contagious" else 0.0
    lam, free_f = {"B2_true": (lam_true, True), "B2_static": (0.0, True), "B2_newton": (0.0, False)}[law]
    p = np.r_[tr["X"][epoch, ic], tr["V"][epoch, ic], np.log(tr["Q"][epoch, ic] / tr["m"][ic]) if free_f else 0.0]
    free = np.arange(7) if free_f else np.arange(6)
    h = np.r_[np.full(3, 1e-7), np.full(3, 1e-6), 1e-5]
    span = LAST - epoch

    def res_batch(Pm):
        lo, la = propagate(tr, epoch, Pm, lam, span)
        return [resid_angles(lo[rel, j], la[rel, j], olo, ola) for j in range(len(Pm))]

    r0 = res_batch(p[None])[0]; cost = (r0 ** 2).sum(); lm = 1e-3
    for it in range(25):
        Pm = np.vstack([p] + [p + h[j] * np.eye(7)[j] for j in free])
        rs = res_batch(Pm); r0 = rs[0]; cost = (r0 ** 2).sum()
        J = np.column_stack([(rs[1 + n] - r0) / h[j] for n, j in enumerate(free)])
        A = J.T @ J; g = J.T @ r0; improved = False
        for _ in range(8):
            dp = np.zeros(7); dp[free] = -np.linalg.solve(A + lm * np.diag(np.diag(A)), g)
            r1 = res_batch((p + dp)[None])[0]; c1 = (r1 ** 2).sum()
            if c1 < cost:
                p, cost, lm, improved = p + dp, c1, lm / 5, True
                break
            lm *= 10
        if not improved or np.abs(dp).max() < 1e-12:
            break
    rms = np.sqrt(cost / len(r0))
    cov = np.zeros((7, 7)); cov[np.ix_(free, free)] = np.linalg.pinv(A) * cost / max(1, len(r0) - len(free))
    sets = np.vstack([p, np.random.default_rng(4).multivariate_normal(p, cov, 20)])
    lo, la = propagate(tr, epoch, sets, lam, END - epoch)
    rows = [(dd, lo_, la_) for dd, b, lo_, la_ in load_truth_pred(world, "predictions_2.csv") if b == c]
    days = np.array([r[0] for r in rows]) - epoch
    E = angsep(lo[days], la[days], np.array([r[1] for r in rows])[:, None], np.array([r[2] for r in rows])[:, None])
    return dict(world=world, label=label, comet=c, law=law, n_obs=len(d), fit_rms_arcmin=float(rms),
                pull_factor=float(np.exp(p[6])) if free_f else 1.0,
                true_pull_factor_at_first_sighting=float(tr["Q"][epoch, ic] / tr["m"][ic]),
                best_fit_median_deg=float(np.median(E[:, 0])), cloud_median_deg=float(np.median(np.median(E[:, 1:], 0))))


def main():
    jobs = [(w, l, c, law) for l, w in KEY.items() for c in COMETS for law in ["B2_true", "B2_static", "B2_newton"]]
    with ProcessPoolExecutor(9) as ex:
        out = list(ex.map(fit_comet, jobs))
    json.dump(out, open(os.path.join(ROOT, "private", "baseline_stage2.json"), "w"), indent=1)
    for l, w in KEY.items():
        print(f"\n{l} ({w}): median error over years 101-120, degrees (best fit / cloud); fit rms in arcmin")
        print(f"  {'comet':5s} {'B2_true':>22s} {'B2_static':>22s} {'B2_newton':>22s}")
        for c in COMETS:
            cells = []
            for law in ["B2_true", "B2_static", "B2_newton"]:
                r = next(x for x in out if x["label"] == l and x["comet"] == c and x["law"] == law)
                cells.append(f"{r['best_fit_median_deg']:7.2f}/{r['cloud_median_deg']:6.2f} ({r['fit_rms_arcmin']:5.1f}')")
            print(f"  {c:5s} " + " ".join(f"{x:>22s}" for x in cells))


if __name__ == "__main__":
    main()
