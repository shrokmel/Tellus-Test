"""Integrator accuracy: run the same system at three tolerances and compare."""
import numpy as np
from contagious import build_system
from world import run
names, m, q, X, V = build_system(60, 1)
R = {}
for rt in [1e-11, 1e-12, 1e-13]:
    for lam in [0, 50]:
        R[(rt, lam)] = run(names, m, q, X, V, lam=lam, years=100, sample_days=10, rtol=rt, tag=f"conv_{rt}_{lam}")
        print(f"rtol {rt:g} lam {lam}: {R[(rt,lam)]['seconds']:.1f}s  {R[(rt,lam)]['log']}")
for lam in [0, 50]:
    for a, b in [(1e-11, 1e-13), (1e-12, 1e-13)]:
        dX = np.linalg.norm(R[(a, lam)]["X"] - R[(b, lam)]["X"], axis=2)
        df = np.abs(R[(a, lam)]["Q"][-1] - R[(b, lam)]["Q"][-1]) / m
        print(f"lam {lam}: rtol {a:g} vs 1e-13: planets max {dX[:,1:7].max():.1e} AU; comets >0.01 AU apart at 100 yr: {(dX[-1,7:]>0.01).sum()}; max pull-factor diff {df[7:].max():.1e}")
