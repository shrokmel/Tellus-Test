"""Pick the 12 comets, build the Tellus system, run both worlds for 100 and 1,000 years."""
import json, numpy as np
from system import assemble, random_comets
from world import run, energy, elements
rows = json.load(open("runs/scan_events.json"))
cands = random_comets(300, seed=7)
with_pass = [53, 268, 70, 3, 167, 279, 217]
quiet = [r for r in rows if not r["events"] and r["rmin_star"] < 2.5]
quiet = sorted(quiet, key=lambda r: r["f0"])
pick_q = [quiet[int(round(x))]["cand"] for x in np.linspace(0, len(quiet) - 1, 5)]
chosen = with_pass + pick_q
print("chosen candidates:", chosen)
names, m, q, X, V = assemble(1, [cands[i] for i in chosen])
json.dump(dict(chosen=chosen, names=names), open("runs/final_choice.json", "w"))
for lam, tag in [(50, "final_contagious"), (0, "final_newton")]:
    r = run(names, m, q, X, V, lam=lam, years=100, sample_days=1, rtol=1e-12, tag=tag)
    np.savez_compressed(f"runs/{tag}.npz", **{k: v for k, v in r.items() if k not in ("log",)})
    E = energy(r)
    print(f"\n{tag}: {r['seconds']:.1f}s ({r['log']}); energy rel change {(E[-1]-E[0])/abs(E[0]):+.1e}")
    f = r["Q"] / m
    for c in range(8, len(m)):
        a0 = elements(r, c, 0)[0]; a1 = elements(r, c, -1)[0]
        print(f"  {names[c]:4s} pull factor {f[0,c]:.3f} -> {f[-1,c]:.3f}   orbit size {a0:6.2f} -> {a1:6.2f} AU")
    am0 = elements(r, 2 + 1, 0, centre=2); am1 = elements(r, 3, -1, centre=2)
    print(f"  moon around Tellus: a {am0[0]:.6f} -> {am1[0]:.6f} AU, e {am0[1]:.4f} -> {am1[1]:.4f}; moon pull factor end {f[-1,3]:.9f}, Tellus {f[-1,2]:.9f}")
