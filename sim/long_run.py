"""Stability: run the final system for 1,000 years in both worlds."""
import json, numpy as np
from system import assemble, random_comets
from world import run, elements
ch = json.load(open("runs/final_choice.json"))
cands = random_comets(300, seed=7)
names, m, q, X, V = assemble(1, [cands[i] for i in ch["chosen"]])
for lam, tag in [(0, "newton"), (50, "contagious")]:
    r = run(names, m, q, X, V, lam=lam, years=1000, sample_days=10, rtol=1e-12, tag=f"final1000_{tag}")
    f = r["Q"] / m
    print(f"== {tag}, 1,000 years: {r['seconds']:.1f}s")
    for i in [1, 2, 4, 5, 6, 7]:
        a0, e0, _ = elements(r, i, 0); a1, e1, _ = elements(r, i, -1)
        print(f"  {names[i]:3s} a {a0:7.4f} -> {a1:7.4f}  e {e0:.4f} -> {e1:.4f}  pull factor {f[-1, i]:.9f}")
    am = [elements(r, 3, k, centre=2) for k in range(0, len(r["t"]), 3650)]
    print("  moon around Tellus every 100 yr, a:", " ".join(f"{x[0]:.5f}" for x in am))
    print("  moon e:", " ".join(f"{x[1]:.3f}" for x in am))
    fates = []
    for c in range(8, len(m)):
        a1 = elements(r, c, -1)[0]
        fates.append(f"{names[c]}:{'escaped' if not np.isfinite(a1) else f'{a1:.1f}'}/{f[-1, c]:.2f}")
    print("  comets (orbit size AU / pull factor at 1,000 yr):", ", ".join(fates))
