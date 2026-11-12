"""predict_2.py -- predictions for all bodies, years 101-120 (fills predictions_2.csv).
Uses only the laws and constants of ANSWER_2.md, section 2:
  Laws 1 and 3 (stage 1, unchanged): frame; M about Tellus                       -> predict.py
  Law 2 (same equation, constants refitted with Tellus' mass m_T != 0): N-body S,Tellus,P1..P5 -> predict.py code
  Law 4 (new): a faint body Ci moves in the field of S, Tellus and P1..P5 with
        d2r/dt2 = kappa_i * sum_j G m_j (r_j - r)/|r_j - r|^3  -  sum_{j != S} G m_j r_j/|r_j|^3
     (heliocentric r; j runs over S [at the origin], Tellus, P1..P5; kappa_i = 1 - beta_i is the body's own
      constant; no other force).  beta_i is piecewise constant: the values below hold after each body's
      last close encounter.
Run:  python predict_2.py
"""
import csv, math
import numpy as np
from scipy.integrate import solve_ivp
import predict as PR          # stage-1 laws and constants (unchanged)

# ---- Law 2 constants, refitted in stage 2 with Tellus' mass (overrides the stage-1 values in predict.py) ----
PR.MU_S = 0.000295846623432512
PR.GM = [8.875403557574778e-10, 1.6870535433300474e-10, 2.7847366215846803e-10, 4.429881178946429e-08, 2.96075243502887e-07, 7.361804764707988e-08]                # GM of T, P1..P5 [L^3/day^2]
PR.X0 = [
    [-0.39427492779104417, -0.9382586701591791, -2.4788743493008835e-05],
    [-0.3125547352535995, -0.30167625499677453, -3.5252971841571567e-06],
    [-0.7147922004209684, 1.63305823887885, 4.2255002974181964e-05],
    [-3.0564423893019383, 1.3649701718651674, 0.05306202883402156],
    [6.460515140610338, -0.028485270938687423, -0.032006542329589924],
    [3.3233313494872707, 11.027411653723954, 0.0001752145109602163]]
PR.V0 = [
    [0.015517484898406611, -0.006687360800559479, -1.7847850177955466e-07],
    [0.017477142339945875, -0.019979125247511095, -1.9403581300883658e-07],
    [-0.011452845264201959, -0.005168155280009651, -9.348256995049019e-08],
    [-0.004117297566362069, -0.008373836406802188, 0.00019967810868422217],
    [0.00021088049390257329, 0.006626428699538542, -0.0001116515449814765],
    [-0.00490180861169432, 0.0012657941030181038, 0.0001766801335408227]]
MU = PR.MU_S
GMP = np.array(PR.GM)
FAINT = {
    'C1': dict(t0=8245.0, state=[-1.1886593861327817, 1.0225476224358698, 0.05556595753238944, -0.02223906287371469, 3.774911035868435e-05, -0.00046330687409200366], beta=-0.01389242716772994),
    'C2': dict(t0=33237.0, state=[-2.237417461197135, -1.0990299667853365, 0.07082076962143544, 0.005797090788055514, -0.010715984234431715, 0.004802961751860533], beta=6.271703758876542e-05),
    'C3': dict(t0=14212.0, state=[1.9143232610815017, -1.0810291394579805, 0.3718878604660159, 0.014214054725223477, 0.008459124469525028, 0.00145076747355194], beta=-0.06346273350050617),
    'C4': dict(t0=35935.0, state=[-0.3875860944629295, 1.5253923951231991, -0.02628201971830218, -0.012884290220295126, 0.00236127779970349, -1.759754665438632e-05], beta=8.229915295148944e-05),
    'C5': dict(t0=34065.0, state=[1.8035697025415187, -1.9198257265456893, -0.07283844442008348, 0.008584846413568601, 0.008139781389789873, 0.00411652934808097], beta=-0.0010126105488754442),
    'C6': dict(t0=34740.0, state=[-0.03204878566920739, -2.759773176194422, -0.6708063242308968, 0.011506344377911133, -0.0005364929650412759, 0.0017001735046026573], beta=0.03515201383551987),
    'C7': dict(t0=36112.0, state=[0.15288569914568031, 1.4121527668428142, -0.2773564499564551, -0.012828245675205205, -0.0013435499273879291, -0.0024511989699116405], beta=-0.0001271336999060671),
    'C8': dict(t0=36123.0, state=[-0.756476016314831, 0.7029822898078553, 0.32045194685758704, -0.011226266611515504, -0.01502976462609772, -0.0004433229691279286], beta=0.29735364934128333),
    'C9': dict(t0=33943.0, state=[-1.0427569876947582, 0.7389181159422937, -0.030535435507154193, -0.007773785222467188, -0.01705697250267892, 0.004623140841808296], beta=0.08181442552922828),
    'C10': dict(t0=27114.0, state=[0.03638771149039295, -0.856755613108223, 0.06173672666810495, 0.025038411581548802, -0.004319873757867915, -0.005607841403389235], beta=-0.03257595035496772),
    'C11': dict(t0=35629.0, state=[-2.220828209502211, 0.08372913142856472, 0.058064433172227564, -0.007599105135812233, -0.012061906299910705, -0.0002459921738597197], beta=-0.1155290572609629),
    'C12': dict(t0=34671.0, state=[-1.743831314572344, -1.355047619175082, -0.045487494543137955, 0.010436810703464, -0.01324851272739797, 0.004750587100370504], beta=-0.2970145148859601),
}             # name -> dict(t0 [day], state[6] at t0 (heliocentric ecliptic, L and L/day), beta)

def planet_table(nd, h=0.25):
    x, v = np.array(PR.X0, float), np.array(PR.V0, float); m = np.array(PR.GM, float)
    N = int(round(nd / h)); X = np.empty((N + 1, 6, 3)); V = np.empty((N + 1, 6, 3)); X[0] = x; V[0] = v
    for s in range(1, N + 1):
        a1 = PR.acc(x, MU, m); x2 = x + 0.5*h*v; v2 = v + 0.5*h*a1
        a2 = PR.acc(x2, MU, m); x3 = x + 0.5*h*v2; v3 = v + 0.5*h*a2
        a3 = PR.acc(x3, MU, m); x4 = x + h*v3; v4 = v + h*a3
        a4 = PR.acc(x4, MU, m)
        x = x + h/6*(v + 2*v2 + 2*v3 + v4); v = v + h/6*(a1 + 2*a2 + 2*a3 + a4)
        X[s] = x; V[s] = v
    return X, V, h

def make_rhs(X, V, H, kappa):
    def planets(t):
        k = min(int(t // H), len(X) - 2); s = t / H - k
        h00 = 2*s**3 - 3*s**2 + 1; h10 = s**3 - 2*s**2 + s; h01 = -2*s**3 + 3*s**2; h11 = s**3 - s**2
        return h00*X[k] + h10*H*V[k] + h01*X[k+1] + h11*H*V[k+1]
    def rhs(t, y):
        r, v = y[:3], y[3:]; rp = planets(t)
        rn = math.sqrt(r @ r)
        a = -MU * r / rn**3
        d = rp - r
        a = a + (GMP[:, None] * d / np.sum(d*d, 1)[:, None]**1.5).sum(0)
        ind = (GMP[:, None] * rp / np.sum(rp*rp, 1)[:, None]**1.5).sum(0)
        return np.r_[v, kappa * a - ind]
    return rhs

def faint_positions(name, X, V, H, days):
    c = FAINT[name]
    rhs = make_rhs(X, V, H, 1.0 - c['beta'])
    te = np.unique(np.asarray(days, float))
    out = {}
    for fwd in (True, False):
        sel = te[te >= c['t0']] if fwd else te[te < c['t0']][::-1]
        if len(sel) == 0: continue
        s = solve_ivp(rhs, (c['t0'], sel[-1]), np.array(c['state'], float), method='DOP853',
                      rtol=1e-11, atol=1e-14, t_eval=sel)
        for t, y in zip(s.t, s.y.T): out[int(round(t))] = y[:3]
    return out

def main():
    rows = list(csv.DictReader(open('predictions_2.csv')))
    days = np.array([int(r['day']) for r in rows]); nd = int(days.max()) + 1
    P = PR.nbody(nd); Mg = PR.moon(nd)
    X, V, H = planet_table(nd)
    names = ['T', 'P1', 'P2', 'P3', 'P4', 'P5']
    fpos = {n: faint_positions(n, X, V, H, days[[r['body'] == n for r in rows]]) for n in FAINT}
    for r, d in zip(rows, days):
        b = r['body']
        if b == 'S': g = -P[d, 0]
        elif b == 'M': g = Mg[d]
        elif b in names: g = P[d, names.index(b)] - P[d, 0]
        else: g = fpos[b][d] - X[int(round(d / H)), 0]
        lo, la = PR.to_obs(g)
        r['lon'] = f'{lo:.5f}'; r['lat'] = f'{la:.5f}'
    with open('predictions_2.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['day', 'body', 'lon', 'lat']); w.writeheader(); w.writerows(rows)
    print('wrote', len(rows), 'predictions')

if __name__ == '__main__':
    main()
