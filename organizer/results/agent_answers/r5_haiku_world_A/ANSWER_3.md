# ANSWER_3

## 1. Laws used for N

Units: angles in radians internally, time in days, distance in AU. The same geometry as stages 1-2.

**Law G (gravitation, Sun).** The Sun's gravitational parameter (stage 1 Law 1) is GM_S = NE^2 = (0.01720207)^2 = 2.95910e-4 AU^3/d^2. Each body's orbit obeys n^2 a^3 = K * GM_S, where K = GM_body/GM_S (K = 1 for stage 1 planets and Earth).

**Law N (new body).** N moves on a heliocentric Kepler orbit with K_N = GM_N / GM_S:
  r_N(t) = Kepler(a_N, e_N, i_N, Omega_N, omega_N, M0_N, n_N; t),  a_N = (K_N * NE^2 / n_N^2)^(1/3).
Its geocentric direction is (r_N - r_E)/|r_N - r_E|, with the Earth orbit of stage 1 (Law 2: e=0.01997902, i=0.37877084, Omega=3.14168505, omega=-1.50264368, M0=15.2238328, n=NE, a=1 AU).

Constants (least-squares fit to all 207 sightings, days 36233-36525; 1 sigma from the covariance matrix; rms residual 0.81 arcmin, the noise floor):

| constant | value | 1 sigma |
|---|---|---|
| e_N | 0.7417993 | 0.0000096 |
| i_N [rad] | 0.5777532 | 0.0000267 |
| Omega_N [rad] | 0.3988447 | 0.0000532 |
| omega_N [rad] | 1.0533917 | 0.0002246 |
| M0_N [rad] (mean anomaly at day 0 of the record) | 2.2066958 | 0.0005075 |
| n_N [rad/d] | 0.00063017 | 0.0000000 (less than 1e-9) |
| K_N = GM_N / GM_S | 1.326294 | 0.000199 |

Derived: a_N = 9.9608 AU; period T_N = 9970.6 d (27.3 yr); perihelion q = 2.572 AU; aphelion Q = 17.350 AU.

**What changed relative to stage 2.** The stage 2 common-GM Kepler fit (K=1) gives 58 arcmin rms for N. A free semi-major-axis fit gives 0.81 arcmin rms, which is the noise floor. The data therefore prefer K_N = 1.326 over K = 1. The change from stage 2 is the one constant: K_N, not the form of the law. The Kepler form itself describes N with no residual structure over the 293-day arc.

## 2. Predictions

`predictions_3.csv` (1461 rows, days 36530 to 43830, the 20 years after the last observation) is filled by `predict_3.py`, which uses only Law G, Law N and the Earth constants above. In-sample check against observations_3.csv: lon rms 0.97 arcmin, lat rms 1.02 arcmin, max 3.2 arcmin.

## 3. N's path

- **Now (day 36525, year 100):** N is at heliocentric distance 2.98 AU, elongation 65 deg from the Sun, at lon 116.6, lat 21.7, moving at 0.47 deg/day on the sky (fast, because it is close to perihelion).
- **Perihelion** was at about day 36381 (year 99.6, a few months before the last sighting), at q = 2.57 AU. The next perihelion is at about day 46350 (year 127), outside the window.
- **Outbound:** N recedes toward aphelion Q = 17.35 AU, at about day 41,400 (year 113). Its distance is 5.4 AU at year 101, 12.6 AU at year 105, 16.7 AU at year 110, and 14.3 AU at year 120.
- **Sky motion** slows sharply: 0.47 deg/day now, 0.09 deg/day at year 101, 0.06 deg/day at year 105, about 0.05 deg/day at year 110-115, and 0.054 deg/day at year 120.
- **Position in the sky:** lon increases from 117 to 287 deg over the window and lat falls from +22 to -36 deg. Elongation from the Sun peaks at about 138 deg near year 105 and falls back to 111 deg by year 120.
- **Possible interaction:** N crossed the radii of P3 (3.3 AU) and P4 (6.2 AU) during years 99.6-101, before the forecast window, so no close pass with them is expected inside it. Its outbound crossing of P5's radius (11.47 AU) comes at about year 104-105, inside the window. A perturbation there is the main risk to the prediction (see section 4).

## 4. Judgement: expected error

The covariance of the fit alone gives about 0.01 deg at 1, 5 and 20 years, but that only measures fit noise, not model error. My estimate, including model error:

- **1 year:** about 0.05 deg. Expected error is small; N is still far from the planets.
- **5 years:** about 0.3 deg. N crosses P5's orbit radius during this span, and a close pass with P5 could shift the direction.
- **20 years:** about 2 deg (range 1 to 5 deg). Errors accumulate through the planetary passes, and the stage 2 faint-body failures show that a single Kepler orbit can drift. The 20-year prediction is the least reliable.

**What would prove me wrong:**
- Any sighting after day 36525 that differs from the predicted lon/lat by more than about 0.1 deg in the first year, or 0.5 deg by year 5.
- A residual that grows steadily over 1-2 years (secular change in the elements, not a one-off jump).
- A systematic offset that depends on the planets' positions (for example, a shift when N passes near P5), which would show a perturbation the stage 1 planets should produce.
- An arc that gives a K_N different from 1.326 when fitted separately: if the first and second halves of the 293-day arc give different K_N, the single-GM Kepler form is wrong for N.

**Confidence:** Kepler form for N over the observed arc, about 80%. K_N = 1.326 as a constant of the orbit, about 60% (it could be absorbed by an unmodelled perturbation during the arc). Predictions to 1 year, about 75%; to 20 years, about 40%.
