# ANSWER_3 - Stage 3: the new body N

## 1. Laws used for N

**Unchanged from stage 2.** N's 207 sightings (days 36233-36525, 292 days, rms residual **1.03 arcmin**, i.e. at the noise) are
fitted by exactly the stage 2 laws, so nothing was revised:

- Law A (massive bodies S, Tellus, M, P1..P5, inverse-square Newtonian gravity, GM_S = 2.9592e-4 AU^3/day^2, frame angle
  EPS = 21.7015 deg): all constants as in ANSWER_1.md / ANSWER_2.md, exactly as coded in `predict_3.py`
  (m_P1 = 2.0e-6, m_P2 = 6.5e-7, m_P3 = 1.4956e-4, m_P4 = 1.0000e-3, m_P5 = 2.498e-4, m_E = 3.0e-6, m_M = 3.7e-8).
- Law B (faint body): N is a massless particle whose acceleration is the Newtonian field of all massive bodies times a
  per-body constant k = 1 - beta:

      d^2 x/dt^2 = (1 - beta_N) * sum_{j in S,E,M,P1..P5} GM_j (r_j - x)/|r_j - x|^3 .

The only things that had to be determined for N are its beta and its six orbital elements. Fitted (200 random starts all
converge to one orbit; osculating heliocentric elements at day 36525 with mu = GM_S (1 - beta); angles in the ecliptic frame; formal 1-sigma):

| constant | value |
|---|---|
| a (AU) | 10.145 +- 0.104 |
| e | 0.74655 +- 0.0026 |
| i (deg) | 53.681 +- 0.002 |
| Omega (deg) | 15.258 +- 0.003 |
| omega (deg) | 70.651 +- 0.017 |
| M0 at day 36525 (deg) | 5.061 +- 0.081 |
| beta (k = 1 - beta) | -0.3221 +- 0.0022 (k = 1.3221) |

(a and e are correlated with beta; the resulting orbital period is ~28 yr, perihelion distance q = 2.57 AU, perihelion about day 36500.)
What is new: nothing in the laws. beta_N = -0.322 is slightly outside the range of the twelve earlier faint bodies (-0.30 to +0.30),
which is not a problem for a law with a free constant per body, but it means k is not confined to 0.70-1.30.

## 2. Predictions
`predict_3.py` (about 50 s) integrates law A from day 0, then N under law B from the elements above at day 36525 to day 43830, and
writes lon, lat (catalogue frame, deg) for the 1461 days of `predictions_3.csv`.

## 3. N's path
- N is a high-inclination (54 deg) comet seen near the end of a perihelion passage: at the last sighting (day 36525) it was at
  catalogue lon 113.7 deg, lat +22.5 deg, 3.0 AU from S and 3.3 AU from Tellus, moving east at ~0.4 deg/day and heading away
  from S and Tellus.
- It keeps moving east and south: lon 119 / lat +21 on day 36530, lon 179 / lat +6 at day 36890 (1 yr), lon 217 / lat -22 at
  day 38350 (5 yr). Its speed on the sky falls from 0.4 deg/day to ~0.05 deg/day around day 38000 as it recedes (distance
  from S 9-10 AU through days 38000-39200).
- Key event: around days 36800-37100 (year ~101) N passes within ~1.4 AU of P4 (mass 1e-3 of S). That strong
  deflection shortens its period from ~28 yr to roughly 12 yr, so it will turn round: it reaches its far point (~10 AU) around
  day 38700, then falls back toward S.
- A new perihelion passage at ~2.2 AU from S around day 41000 (year ~112). Between about days 40200 and 41300 N speeds up again
  (0.2-0.4 deg/day), swings far south (lat down to -59 deg near day 40600) and sweeps across lon ~250 -> 359 -> 170 deg. After that it
  recedes again: lon 205 / lat -9 at day 41975, lon 233 / lat -34 at day 43825.
- It should be faint and unseen for long stretches (only recorded near perihelion), probably visible again from about day 40300 to 41700.

## 4. Self-assessment
Expected error of the predicted direction (angle on the sky, from sampling the fit covariance through the full N-body integration; median, with the
95th percentile in brackets):

| after | error |
|---|---|
| 1 year (day ~36890) | 0.06 deg (0.14) |
| 5 years (day ~38350) | 2.8 deg (9) |
| 20 years (day ~43825) | 13 deg (34) |

(and ~16 deg median, 56 deg 95th percentile, during the fast phase near day 40200.) The growth is because the close pass of P4
amplifies small orbit errors, and because period/timing errors move the second perihelion swing.
These errors assume law B (including the planets' pull scaled by the same k, which held to 0.4% for three other bodies) is right. If a
close encounter with P4 is actually closer or farther than computed, the later path (after year ~102) can be completely different.

What would prove me wrong: N being more than ~0.3 deg from the prediction by day ~37000 (it would mean the P4 pass or k is different);
N not turning back toward S by day ~39000 (would mean no strong deflection); a change in N's apparent direction larger
than a few degrees from parallax-free motion at large distance, or any sign that k for N changes in time.
