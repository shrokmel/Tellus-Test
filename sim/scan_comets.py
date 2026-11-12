"""Scan candidate comets (as test particles) in the contagious world; record their close passes."""
import json, numpy as np
from system import assemble, random_comets
from world import run
cands = random_comets(300, seed=7)
names, m, q, X, V = assemble(1, cands)
tp = np.array([0] * 8 + [1] * len(cands))
r = run(names, m, q, X, V, lam=50, d=0.1, years=100, sample_days=5, rtol=1e-12, tp=tp, tag="scan")
print(f"scan: {r['seconds']:.1f}s {r['log']}")
f = r["Q"] / m; t = r["t"]
rows = []
for c in range(8, len(m)):
    i = c - 8
    df = np.abs(np.diff(f[:, c]))
    # events: samples where the pull factor moves by more than 0.002 in 5 days
    ev = []
    idx = np.where(df > 0.002)[0]
    groups = np.split(idx, np.where(np.diff(idx) > 4)[0] + 1) if len(idx) else []
    for g in groups:
        k0, k1 = g[0], g[-1] + 1
        dist = np.linalg.norm(r["X"][k0:k1 + 1, :8] - r["X"][k0:k1 + 1, c][:, None], axis=2)
        j = int(np.argmin(dist.min(0)))
        ev.append(dict(year=round(float(t[k0]), 2), body=names[j], dmin=round(float(dist.min()), 4),
                       f_before=round(float(f[k0, c]), 4), f_after=round(float(f[k1, c]), 4)))
    rmin_star = float(np.linalg.norm(r["X"][:, c] - r["X"][:, 0], axis=1).min())
    rows.append(dict(cand=i, f0=round(float(f[0, c]), 4), f_end=round(float(f[-1, c]), 4), events=ev,
                     rmin_star=round(rmin_star, 3)))
json.dump(rows, open("runs/scan_events.json", "w"), indent=1)
n_ev = sum(1 for r_ in rows if r_["events"])
print(f"candidates with at least one charge-changing pass: {n_ev} of {len(rows)}")
big = [r_ for r_ in rows if any(abs(e['f_after'] - e['f_before']) > 0.05 for e in r_['events'])]
print(f"with a pass changing the pull factor by > 0.05: {len(big)}")
for r_ in big[:40]:
    print(r_["cand"], r_["f0"], [(e["year"], e["body"], e["dmin"], e["f_before"], e["f_after"]) for e in r_["events"]])
