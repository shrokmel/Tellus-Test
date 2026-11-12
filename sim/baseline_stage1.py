"""Stage 1 baselines: predict S, P1-P5 and M in years 101-120 from the 100 years of bright-body sightings.

  B1_newton : the true law for these bodies (Newtonian gravity, all pulling on each other). Fitted parameters: the
              state of P1, Tellus, M, P2-P5 relative to S at day 0 (42) and their log-masses (7); S's mass sets the
              units. Idealised: the search starts at the truth.
  B1_kepler : the same, but the planets do not pull on each other (each on its own ellipse around S; M feels S
              and Tellus). Shows how much the planets' mutual pulls matter. No masses fitted (42 parameters).
Fit: Gauss-Newton / Levenberg-Marquardt on angular residuals; Jacobian by finite differences, one integration per
parameter, run in parallel. Uncertainty: 20 parameter sets drawn from the fit's covariance (B1_newton only).
Score: great-circle error against the true direction on the prediction days (every 30 days, years 101-120).
Writes ../private/baseline_stage1.json.
"""
import json, os, sys, numpy as np
from concurrent.futures import ProcessPoolExecutor
from baseline_lib import ROOT, LAST, END, KEY, truth_run, load_obs, load_truth_pred, integrate, resid_angles, angsep, directions_from

BODIES = ["S", "P1", "T", "M", "P2", "P3", "P4", "P5"]
FIT = [b for b in BODIES if b != "S"]
OBS_BODIES = ["S", "P1", "P2", "P3", "P4", "P5", "M"]
G_WORKERS = 8
VARIANTS = sys.argv[1:] or ["newton", "kepler"]


def setup(world):
    tr = truth_run(world)
    idx = [tr["names"].index(b) for b in BODIES]
    return dict(X=tr["X"][0, idx], V=tr["V"][0, idx], m=tr["m"][idx], Q=tr["Q"][0, idx])


def build(base, p, variant):
    """Parameters -> integrator input. p = 6 per FIT body (relative to S) [+ 7 log-masses for newton]."""
    X = base["X"].copy(); V = base["V"].copy(); m = base["m"].copy()
    for k, b in enumerate(FIT):
        i = BODIES.index(b)
        X[i] = X[0] + p[6 * k:6 * k + 3]; V[i] = V[0] + p[6 * k + 3:6 * k + 6]
    if variant == "newton":
        for k, b in enumerate(FIT):
            m[BODIES.index(b)] = np.exp(p[42 + k])
    tp = np.zeros(len(BODIES), int)
    if variant == "kepler":
        for b in ["P1", "P2", "P3", "P4", "P5", "M"]:
            tp[BODIES.index(b)] = 1
    return X, V, m, m.copy(), tp            # pull factor 1: q = m


def simulate(args):
    base, p, variant, days = args
    X, V, m, q, tp = build(base, p, variant)
    r = integrate(BODIES, m, q, X, V, 0.0, days, tp=tp)
    iT = BODIES.index("T")
    return {b: directions_from(r, iT, BODIES.index(b)) for b in OBS_BODIES}


def residuals(dirs, obs, upto=LAST):
    out = []
    for b in OBS_BODIES:
        d, olo, ola = obs[b]
        keep = d <= upto; d, olo, ola = d[keep], olo[keep], ola[keep]
        out.append(resid_angles(dirs[b][0][d], dirs[b][1][d], olo, ola))
    return np.concatenate(out)


def p_truth(base, variant):
    p = []
    for b in FIT:
        i = BODIES.index(b)
        p += list(base["X"][i] - base["X"][0]) + list(base["V"][i] - base["V"][0])
    if variant == "newton":
        p += [np.log(base["m"][BODIES.index(b)]) for b in FIT]
    return np.array(p)


def steps(n):
    s = np.r_[np.tile(np.r_[np.full(3, 1e-7), np.full(3, 1e-6)], 7)]
    return np.r_[s, np.full(n - len(s), 1e-4)] if n > len(s) else s


def fit(base, obs, variant, ex, iters=8, upto=LAST, p0=None):
    p = p_truth(base, variant) if p0 is None else p0.copy(); h = steps(len(p)); lam_lm = 1e-3
    r0 = residuals(simulate((base, p, variant, upto)), obs, upto); cost = (r0 ** 2).sum()
    for it in range(iters):
        cols = list(ex.map(simulate, [(base, p + h[j] * np.eye(len(p))[j], variant, upto) for j in range(len(p))]))
        J = np.column_stack([(residuals(c, obs, upto) - r0) / h[j] for j, c in enumerate(cols)])
        A = J.T @ J; g = J.T @ r0
        improved = False
        for _ in range(8):
            dp = -np.linalg.solve(A + lam_lm * np.diag(np.diag(A)), g)
            r1 = residuals(simulate((base, p + dp, variant, upto)), obs, upto); c1 = (r1 ** 2).sum()
            if c1 < cost:
                p, r0, cost, lam_lm, improved = p + dp, r1, c1, lam_lm / 5, True
                break
            lam_lm *= 10
        rms = np.sqrt(cost / len(r0))
        print(f"    {variant} upto day {upto} iter {it}: rms {rms:.4f}'  step {np.abs(dp).max():.1e}", flush=True)
        if not improved or np.abs(dp).max() < 1e-12:
            break
    cov = np.linalg.pinv(A) * cost / (len(r0) - len(p))
    return p, cov, np.sqrt(cost / len(r0))


def score(world, dirs_list):
    rows = load_truth_pred(world, "predictions.csv")
    errs = {b: [] for b in OBS_BODIES}
    for dirs in dirs_list:
        e = {b: [] for b in OBS_BODIES}
        for d, b, lo, la in rows:
            e[b].append(angsep(dirs[b][0][d], dirs[b][1][d], lo, la))
        for b in OBS_BODIES:
            errs[b].append(np.array(e[b]))
    out = {}
    for b in OBS_BODIES:
        E = np.array(errs[b])                       # (sets, days)
        out[b] = dict(best_fit_median_deg=float(np.median(E[0])), best_fit_max_deg=float(E[0].max()),
                      cloud_median_deg=float(np.median(np.median(E[1:], axis=1))) if len(E) > 1 else None)
    return out


def main():
    res = {}
    with ProcessPoolExecutor(G_WORKERS) as ex:
        for label, world in KEY.items():
            base = setup(world); obs = load_obs(label, 1)
            print(f"{label} ({world}): {sum(len(v[0]) for v in obs.values())} sightings", flush=True)
            res[label] = dict(world=world)
            for variant in VARIANTS:
                if variant == "kepler":
                    p = None
                    for upto in [3653, 10958, 21915, LAST]:          # 10, 30, 60, 100 years
                        p, cov, rms = fit(base, obs, variant, ex, iters=12, upto=upto, p0=p)
                else:
                    p, cov, rms = fit(base, obs, variant, ex)
                sets = [p]
                if variant == "newton":
                    sets += list(np.random.default_rng(3).multivariate_normal(p, cov, 20))
                dirs = list(ex.map(simulate, [(base, s, variant, END) for s in sets]))
                sc = score(world, dirs)
                res[label]["B1_" + variant] = dict(fit_rms_arcmin=rms, per_body=sc)
                if variant == "newton":
                    res[label]["B1_newton"]["fitted_masses"] = {b: float(np.exp(p[42 + k])) for k, b in enumerate(FIT)}
                print(f"  B1_{variant}: fit rms {rms:.3f}'; median error years 101-120 (deg): " +
                      ", ".join(f"{b} {v['best_fit_median_deg']:.4f}" for b, v in sc.items()), flush=True)
    path = os.path.join(ROOT, "private", "baseline_stage1.json")
    old = json.load(open(path)) if os.path.exists(path) else {}
    for l in res:
        old.setdefault(l, {}).update(res[l])
    json.dump(old, open(path, "w"), indent=1)


if __name__ == "__main__":
    main()
