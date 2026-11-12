# ANSWER_2 - Stage 2: the faint bodies C1-C12

## 1. Test of the stage 1 laws against all bodies

**Where they hold.**
- Stage 1 bodies (S, M, P1-P5): the stage 1 law (Newtonian inverse-square gravity between S, Tellus, M, P1-P5) still reproduces all
  156,000 old sightings at the noise level: rms 0.99-1.02 arcmin in each angle for every body (worst: M longitude 1.01).
  Only one constant had to change (see below): the mass of P2 (see section 2) and, to compensate, tiny re-adjustments
  of P2's and P3's initial elements (shifts of ~1e-5 AU / 0.002 deg).
- The faint bodies are massless test particles: nothing in the planets' motion needs them (P3-P5 and S stay at noise
  with C1-C12 ignored), and the *form* of the law (inverse-square attraction towards S and all massive bodies) is
  right for them too.
- Frame, EPS = 21.7015 deg, and Tellus's orbit: unchanged.

**Where they fail.** Applied literally (a faint body is a test particle with the same acceleration as a planet) the stage 1 laws
fail for every faint body, by a lot:
- C2, 260-day arc: best possible orbit around S under gravity of S+planets leaves rms 4.7 arcmin (smooth +-17' residuals),
  on a 490-day arc 11.3 arcmin (noise is 1.0). Over several apparitions the best fit is much worse (rms 4 to 20 arcmin).
- C7 (arc of 615 d): best fit rms 100 arcmin; C1: > 2000 arcmin with 2 apparitions.
- In every case the residual is cured by one extra constant per body: the pull of S is weaker or stronger than for the
  planets by a fixed factor. With it the same bodies are fitted at 1.00-1.13 arcmin rms over the *whole* 100 years
  (e.g. C7: 19,086 sightings, rms 1.07; C2: 1.13; C3: 1.03).
- Not needed / ruled out: transverse or out-of-plane extra accelerations (fitted for C2: < 1e-3 of the solar pull);
  a different exponent of r (all fits at noise with exponent 2).

## 2. Laws as equations

**Law A (unchanged from stage 1) - massive bodies.** For i in {S, E (Tellus), M, P1..P5}:

    d^2 r_i/dt^2 = sum_{j != i} GM_j (r_j - r_i) / |r_j - r_i|^3 ,    GM_j = GM_S m_j ,    GM_S = 2.9592e-4 AU^3/day^2.

Initial elements, masses, EPS = 21.7015 deg and the observation equations are exactly those of ANSWER_1.md
(with m_P2 = 6.5e-7 instead of 3e-7, and P2 elements a = 1.699994, e = 0.050003, varpi = 281.1828, L0 = 114.9177;
P3: a = 3.300006, e = 0.039993, varpi = 42.5449, L0 = 151.6882, i = 1.4990, Omega = 118.6446; all as listed in
`predict_2.py`; m_E = 3.0e-6, m_M = 3.7e-8, m_P1 = 2.0e-6 are not measured by the directions, m_P2 = (6.5 +- 1.5)e-7 is
measured weakly by close passes of C3, C7 and others).

**Law B (new) - faint bodies.** Each faint body C_k is a massless particle whose acceleration is the Newtonian gravitational
acceleration multiplied by its own constant factor k_k = 1 - beta_k:

    d^2 x/dt^2 = (1 - beta_k) * sum_{j in S,E,M,P1..P5} GM_j (r_j - x) / |r_j - x|^3 .

(Equivalently: the Sun's pull is reduced to (1-beta_k) of the planet value, *and* so is every planet's pull on the body. A test
with an independent factor kappa on the planetary pull gave kappa/(1-beta) = 1.000-1.004 for C2, C3, C7.) The factor
is a property of each body, between 0.70 and 1.30, not the same for all (beta from -0.30 to +0.30).
The observed direction is that of (x - r_E) exactly as for the other bodies.

**Constants (day 36525 osculating heliocentric elements with mu = GM_S (1-beta); a in AU, angles in deg, M0 mean anomaly, ecliptic frame).**
Uncertainties are formal 1-sigma from the fit covariance (the real ones are several times larger: parameters are correlated and
P1/P2 masses are not exactly known); the last columns are rms (arcmin, both angles, all data) and the number of sightings.

| body | a | e | i | Omega | omega | M0 | beta (k = 1-beta) | k | rms | n |
|---|---|---|---|---|---|---|---|---|---|---|
| C1 | 13.1921 +- 0.0002 | 0.94708 +- 0.00001 | 4.608 +- 0.001 | 343.697 +- 0.004 | 79.863 +- 0.005 | 291.951 | -0.2649 +- 0.0001 | 1.2649 | 0.99 | 449 |
| C2 | 7.4695 +- 0.0002 | 0.77989 +- 0.00001 | 0.264 +- 0.001 | 193.0 +- 0.3 | 8.6 +- 0.3 | 316.933 | 0.0992 +- 0.0001 | 0.9008 | 1.13 | 1176 |
| C3 | 5.7472 +- 0.0001 | 0.69747 +- 0.00001 | 9.622 +- 0.001 | 239.084 +- 0.004 | 36.128 +- 0.003 | 357.732 | -0.2310 +- 0.0001 | 1.2310 | 1.03 | 2222 |
| C4 | 16.2948 +- 0.0003 | 0.91704 +- 0.00001 | 5.111 +- 0.014 | 143.221 +- 0.014 | 223.908 +- 0.015 | 302.247 | 0.2736 +- 0.0001 | 0.7264 | 1.00 | 501 |
| C5 | 4.3702 +- 0.0001 | 0.38735 +- 0.00001 | 18.164 +- 0.001 | 315.409 +- 0.001 | 6.243 +- 0.001 | 77.375 | -0.1721 +- 0.0001 | 1.1721 | 1.12 | 3431 |
| C6 | 5.7432 +- 0.0004 | 0.50018 +- 0.00002 | 15.948 +- 0.002 | 326.439 +- 0.006 | 322.235 +- 0.005 | 351.119 | 0.1441 +- 0.0002 | 0.8559 | 1.00 | 646 |
| C7 | 2.0293 +- 0.0001 | 0.50680 +- 0.00001 | 17.461 +- 0.001 | 223.456 +- 0.001 | 359.039 +- 0.001 | 179.151 | 0.2171 +- 0.0001 | 0.7829 | 1.07 | 19086 |
| C8 | 6.4680 +- 0.0001 | 0.83352 +- 0.00001 | 17.192 +- 0.001 | 48.888 +- 0.002 | 101.910 +- 0.002 | 10.538 | 0.2980 +- 0.0001 | 0.7020 | 1.11 | 1876 |
| C9 | 5.1906 +- 0.0001 | 0.76518 +- 0.00001 | 13.940 +- 0.001 | 150.193 +- 0.002 | 19.306 +- 0.002 | 203.759 | 0.0824 +- 0.0001 | 0.9176 | 1.04 | 2502 |
| C10 | 9.1095 +- 0.0006 | 0.90884 +- 0.00001 | 14.148 +- 0.001 | 108.680 +- 0.005 | 139.662 +- 0.004 | 350.053 | -0.0355 +- 0.0002 | 1.0355 | 1.05 | 762 |
| C11 | 3.2139 +- 0.0001 | 0.54957 +- 0.00001 | 2.576 +- 0.001 | 32.978 +- 0.008 | 57.835 +- 0.008 | 283.436 | -0.1574 +- 0.0001 | 1.1574 | 1.12 | 6350 |
| C12 | 9.4899 +- 0.0005 | 0.76774 +- 0.00001 | 15.780 +- 0.001 | 221.966 +- 0.004 | 357.217 +- 0.005 | 71.049 | -0.2974 +- 0.0002 | 1.2974 | 1.02 | 885 |

(Full-precision values: the COMET table in `predict_2.py`. a is the osculating value at day 36525, which for far-out bodies such
as C4 differs from the long-term mean because the planets perturb it.)
Orbital periods (osculating, mu = GM_S k): C1 ~43 yr, C2 21 yr, C3 12.3 yr, C4 32-77 yr (osculating value varies strongly with epoch), C5 8.5 yr, C6 14.7 yr, C7 3.3 yr, C8 19.5 yr,
C9 12.4 yr, C10 27 yr, C11 5.4 yr, C12 ~25 yr. Perihelion distances 0.7-2.8 AU; the bodies are recorded only
near perihelion (when they are inside ~2-3 AU of S / bright).

Same as stage 1? The planets' law and constants: yes, same (one planet mass, P2, moved by 3e-7 -> 6.5e-7: forced by the faint bodies'
close passes; stage 1 data cannot tell, they are fitted equally well with either after a tiny re-adjustment of the P2 and P3
elements). The new law B is forced by the faint bodies: without the factor k_k no gravitational orbit around S fits even one 490-day
arc of C2 (4.7-11 arcmin vs 1.0 noise), and with it all 12 bodies fit over 100 years.

## 3. Predictions

`predict_2.py` integrates law A from the day-0 state (about 45 s), then law B for the 12 faint bodies from the day 36525 elements
above up to day 43815 (RK4, 1/8-day steps, with the massive bodies' positions from law A), and writes `predictions_2.csv`
(lon, lat in the catalogue frame, degrees) for all 7 + 12 bodies, including the faint bodies on days when they will not be
visible. It does not read the observation files; its constants are exactly those in section 2. (Cross-check: it agrees to
<0.02 arcmin with my C-coded fitting pipeline.)
How well will it do? Hold-out test: refitting only the first 80 years and predicting the last 20 years of actual sightings:
median error 1.2-2.3 arcmin (max 2.8-5.5) for C5, C6, C7, C8, C9, C11, C12; C3: median 2.3, max 14.5; C2: median 27 arcmin
(its fit is the least secure: it has close passes of P2, about 0.25 AU). C1 and C4 have a single apparition in the
prediction window with only 1-2 observed apparitions behind them, so their errors could be tens of arcmin, and
timing errors for C1's perihelion (~day 39700) of a few days are possible.

## 4. Self-assessment

- Law A (stage 1 gravity for massive bodies): 97%. (Unchanged and still at the noise level.)
- Law B in form "acceleration = k_k x Newtonian field, k_k constant per body": 93%. It fits 12 bodies, 40,000 sightings, 100 years, to
  1.0-1.13 arcmin. The numbers k_k have no pattern I can find (range 0.70-1.30), so I do not know why they differ.
- That the planets' pull on a faint body is scaled by the same k (rather than only S's pull): 85% (fitted independent factor
  1.000-1.004 of k for three bodies; for C2 the fit improved from 1.24 to 0.99 with the freed factor, so a small deviation or a
  mass error is possible).
- That faint bodies do not act back on the planets: 90% (not detectable).
- The masses m_P1 = 2e-6 and m_P2 = 6.5e-7: P2 about 70% right within a factor 1.5; P1, Tellus and M masses are unmeasured: 20%.
- Residual excess over noise (rms 1.11-1.13 for C2, C5, C8, C11 against 1.00) I attribute to these mass uncertainties; 60% confident.

What would prove me wrong:
- A faint body at the predicted place at the wrong time or place by more than ~10-30 arcmin when it is next seen (C5 and C11
  are the most secure; C2, C1, C4 least).
- A faint body whose k changes in time (e.g. k drifting between apparitions) or depends on distance from S (a fit with one
  constant k would then fail at large or small r; none seen between r = 0.7 and 25 AU).
- A faint body that deflects a planet at a close pass.
- A well-measured close encounter that needs a planetary pull on the body different from k x the planet's pull.
