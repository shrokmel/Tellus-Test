"""Run the N candidates from design_n.py forward in the full contagious system and report each one's pass."""
import json, numpy as np
from system import assemble, random_comets
from world import run
from observe import visible
ch = json.load(open("runs/final_choice.json")); cands = random_comets(300, seed=7)
D = json.load(open("runs/n_design.json"))
names, m, q, X, V = assemble(1, [cands[i] for i in ch["chosen"]])
names = list(names) + [f"n{i}" for i in range(len(D))]
m = np.r_[m, np.full(len(D), 1e-16)]; q = np.r_[q, np.array([d["f"] for d in D]) * 1e-16]
X = np.vstack([X, [d["x0"] for d in D]]); V = np.vstack([V, [d["v0"] for d in D]])
tp = np.r_[np.zeros(len(names) - len(D), int), np.ones(len(D), int)]
r = run(names, m, q, X, V, lam=50, years=121, sample_days=1, rtol=1e-12, tp=tp, tag="ncheck")
t = r["t"]; f = r["Q"] / m; win = (t >= 36465 / 365.25) & (t <= 100 + 1e-9)
print(f"{r['seconds']:.1f}s")
for i, d in enumerate(D):
    c = 20 + i
    k = np.argmax(np.abs(np.diff(f[:, c]))) + 1
    dist = np.linalg.norm(r["X"][k, 1:8] - r["X"][k, c], axis=1)
    pre = np.abs(np.diff(f[t < 99.8, c])).sum()
    print(f"n{i}: planned pass {d['planet']} yr {d['tp']:.2f} | actual biggest change yr {t[k]:.2f} near {names[1+np.argmin(dist)]} "
          f"({dist.min():.3f} AU); pull factor {f[0,c]:.3f} -> {f[-1,c]:.3f}; change before window {pre:.1e}; "
          f"visible nights in window (no weather) {visible(r['X'][win], 2, 0, c, True).sum()}/60")
