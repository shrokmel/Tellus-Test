"""Design N by working backwards from a chosen close pass.

1. Planets + moon run forward (Newton; planets are identical in both worlds) to get their state at the pass time tp.
2. Candidates: N placed 0.02-0.08 AU from a planet at tp, with a random relative velocity, pull factor in 0.65-0.8
   or 1.2-1.35, on a bound orbit.
3. Run backwards (all velocities reversed, plain Newton: exact, because N meets nothing before the pass) to t = 0.
4. Keep candidates that (a) come no closer than 0.3 AU to the star, any planet or the moon before tp - 0.1 yr,
   (b) are visible by the faint-body rule on >= 8 of the 12 five-day samples in the 60-night window,
   (c) have tp between 100.3 and 104 years.
Output: runs/n_design.json (candidates with their state at t = 0).
"""
import json, numpy as np
from concurrent.futures import ProcessPoolExecutor
from system import planets_and_moon, G
from world import run
from observe import visible
DT = 5 / 365.25
WIN0, WIN1 = 36465 / 365.25, 36525 / 365.25

def planets_run():
    names, m, f, X, V = planets_and_moon(1)
    X = X - (m[:, None] * X).sum(0) / m.sum(); V = V - (m[:, None] * V).sum(0) / m.sum()
    return run(names, m, f * m, X, V, lam=0, years=105, sample_days=5, rtol=1e-12, tag="planets_only")

def batch(b):
    P = dict(np.load("runs/planets_only.npz", allow_pickle=True)); P["names"] = list(P["names"]); rng = np.random.default_rng(500 + b)
    names0, m0 = P["names"], P["m"]
    ks = np.arange(int(np.ceil(100.3 / DT)), int(104 / DT))
    out = []
    for kk in rng.choice(ks, 5):
        tp = P["t"][kk]; Xp, Vp = P["X"][kk], P["V"][kk]
        cand = []
        while len(cand) < 120:
            j = rng.choice([1, 2, 4, 5, 6])          # P1, Tellus, P2, P3, P4
            off = rng.normal(size=3); off *= rng.uniform(0.02, 0.08) / np.linalg.norm(off)
            vrel = rng.normal(size=3); vrel *= rng.uniform(2, 9) / np.linalg.norm(vrel)
            fk = rng.choice([rng.uniform(0.65, 0.8), rng.uniform(1.2, 1.35)])
            x = Xp[j] + off; v = Vp[j] + vrel
            r = np.linalg.norm(x - Xp[0]); w = v - Vp[0]
            if 0.5 * w @ w - G * fk / r >= -0.05: continue      # keep only clearly bound orbits
            cand.append((j, fk, x, v))
        names = list(names0) + [f"n{i}" for i in range(len(cand))]
        m = np.r_[m0, np.full(len(cand), 1e-16)]
        q = np.r_[m0, np.array([c[1] for c in cand]) * 1e-16]
        X = np.vstack([Xp, [c[2] for c in cand]]); V = -np.vstack([Vp, [c[3] for c in cand]])
        tp_flag = np.r_[np.zeros(len(m0), int), np.ones(len(cand), int)]
        R = run(names, m, q, X, V, lam=0, years=tp, sample_days=5, rtol=1e-12, tp=tp_flag, tag=f"nback{b}")
        tf = tp - R["t"]                                   # forward time of each backward sample
        win = (tf >= WIN0 - 1e-9) & (tf <= WIN1 + 1e-9)
        early = tf < tp - 0.1
        for i, (j, fk, x, v) in enumerate(cand):
            c = len(m0) + i
            dmin = np.linalg.norm(R["X"][early][:, :len(m0)] - R["X"][early][:, c][:, None], axis=2).min()
            if dmin < 0.3: continue
            vis = visible(R["X"][win], 2, 0, c, faint=True).sum()
            if vis < 8: continue
            out.append(dict(tp=float(tp), planet=names0[j], f=float(fk), dmin_before=float(dmin), vis=int(vis),
                            x0=R["X"][-1, c].tolist(), v0=(-R["V"][-1, c]).tolist(), t_end_back=float(R["t"][-1])))
    return out

if __name__ == "__main__":
    P = planets_run(); np.savez("runs/planets_only.npz", **{k: v for k, v in P.items() if k != "log"})
    with ProcessPoolExecutor(4) as ex:
        res = [g for gs in ex.map(batch, range(4)) for g in gs]
    print("candidates passing:", len(res), "of 2,400")
    for g in res[:30]:
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in g.items() if k not in ("x0", "v0")})
    json.dump(res, open("runs/n_design.json", "w"), indent=1)
