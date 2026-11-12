"""Build star systems, run them with the C integrator (sky.c), and load the results.

    from world import run
    r = run(names, m, q, X, V, lam=50, d=0.1, years=100, sample_days=1, tp=None, tag="contagious")
    r["t"], r["X"] (samples x N x 3), r["V"], r["Q"], r["m"], r["names"]
"""
import os, subprocess, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
G = 4 * np.pi ** 2


def build_sky(force=False):
    exe = os.path.join(HERE, "sky")
    src = os.path.join(HERE, "sky.c")
    if force or not os.path.exists(exe) or os.path.getmtime(exe) < os.path.getmtime(src):
        subprocess.run(["cc", "-O3", "-march=native", "-ffast-math", "-o", exe, src, "-lm"], check=True)
    return exe


def run(names, m, q, X, V, lam, d=0.1, years=100.0, sample_days=1.0, rtol=1e-10, tp=None, tag="run",
        outdir=os.path.join(HERE, "runs")):
    exe = build_sky()
    os.makedirs(outdir, exist_ok=True)
    N = len(m)
    tp = np.zeros(N, int) if tp is None else np.asarray(tp, int)
    inp = os.path.join(outdir, f"{tag}.in.txt"); out = os.path.join(outdir, f"{tag}.bin")
    with open(inp, "w") as fh:
        fh.write(f"{N} {float(G)!r} {float(lam)!r} {float(d)!r} {float(years)!r} {float(sample_days)!r} {float(rtol)!r}\n")
        for i in range(N):
            fh.write(" ".join([names[i]] + [repr(float(v)) for v in
                     (m[i], q[i], *X[i], *V[i])] + [str(int(tp[i]))]) + "\n")
    t0 = time.time()
    p = subprocess.run([exe, inp, out], capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(p.stderr)
    raw = np.fromfile(out, dtype=np.float64).reshape(-1, 1 + 7 * N)
    st = raw[:, 1:].reshape(-1, 7 * N)
    res = dict(t=raw[:, 0], X=st[:, :3 * N].reshape(-1, N, 3), V=st[:, 3 * N:6 * N].reshape(-1, N, 3),
               Q=st[:, 6 * N:], m=np.asarray(m, float), names=list(names), lam=lam, d=d, tp=tp,
               seconds=time.time() - t0, log=p.stderr.strip())
    return res


def energy(r):
    """Total energy at every sample (gravity between non-test bodies, plus test bodies' energy in the field)."""
    m, Q, X, V, tp = r["m"], r["Q"], r["X"], r["V"], r["tp"]
    N = len(m); iu = np.triu_indices(N, 1)
    keep = ~(tp[iu[0]].astype(bool) & tp[iu[1]].astype(bool))
    E = np.empty(len(r["t"]))
    for k in range(len(r["t"])):
        kin = 0.5 * (m * (V[k] ** 2).sum(1)).sum()
        D = X[k][None] - X[k][:, None]; rr = np.sqrt((D ** 2).sum(-1))[iu]
        E[k] = kin - G * ((Q[k][:, None] * Q[k][None, :])[iu][keep] / rr[keep]).sum()
    return E


def elements(r, i, k, centre=0):
    """Osculating a, e, inclination (deg) of body i around body `centre` at sample k, using the true pull."""
    x = r["X"][k, i] - r["X"][k, centre]; v = r["V"][k, i] - r["V"][k, centre]
    m, Q = r["m"], r["Q"][k]
    mu = G * Q[centre] * Q[i] / m[i] * (1 + m[i] / m[centre])
    R = np.linalg.norm(x); E = 0.5 * v @ v - mu / R
    if E >= 0:
        return np.inf, np.nan, np.nan
    a = -mu / (2 * E); h = np.cross(x, v); e = np.linalg.norm(np.cross(v, h) / mu - x / R)
    return a, e, np.degrees(np.arccos(h[2] / np.linalg.norm(h)))
