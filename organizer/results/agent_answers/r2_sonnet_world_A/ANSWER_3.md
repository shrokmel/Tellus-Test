# Stage 3: the new body N

`python predict_3.py` (needs `predict.py`, `predict_2.py` in this folder) writes `predictions_3.csv` (1461 days, 36530-43830 step 5).

## 1. Laws used for N (unchanged from stage 2)

N's 207 sightings (days 36233-36519) are fitted to noise (median vector residual 1.22 arcmin, same as pure noise) by the stage-2 law for faint bodies, so I did not change the laws:

    d2r/dt2 = -(1 - beta) G m_S r/|r|^3 + sum_{j in {T,P2,P3,P4,P5}} G m_j [ (r_j - r)/|r_j - r|^3 - r_j/|r_j|^3 ]

r is heliocentric in the ecliptic frame (tilt eps = 21.702 deg from the observation frame, as in stage 1), G m_S = n_T^2 = 2.959e-4 AU^3/day^2 (n_T = 0.017202503 rad/day), m_T = 3.0e-6 m_S, m_P2 = 0, m_P3 = 1.5146e-4, m_P4 = 1.0003e-3, m_P5 = 2.491e-4 (times m_S), planets' positions from the stage-1 N-body model (constants in ANSWER_1.md / `predict.py`).

What the data of N determine (fit of the osculating state at epoch day 36379 and beta; 1-sigma from the covariance of the fit):

| constant | value |
|---|---|
| beta_N | -0.3233 +- 0.0017 |
| x, y, z (AU) | 0.470338 +- 9e-5, 1.612196 +- 1.4e-4, 1.947619 +- 2.4e-4 |
| vx, vy, vz (AU/day) | -0.01566352 +- 7.6e-7, -0.00089765 +- 1.1e-6, 0.00443107 +- 1.1e-6 |

(ecliptic frame, epoch 36379 = close to perihelion; a = 10.09 AU, e = 0.745, q = 2.57 AU, i = 53.7 deg, osculating with mu = (1-beta) G m_S.) beta is sharply determined: with beta = -0.30 the fit worsens to 1.31 arcmin, at -0.20 to 3.7, at 0 to 9.7 (so the "extra attraction" seen for several faint bodies in stage 2 is again needed: N's effective solar attraction is 1.32 times Newton's). The solution was unique over a wide grid of beta (-0.6 to 0.9). Nothing in N's sightings changes my mind about the laws; the main uncertainties are the ones inherited from stage 2 (beta being constant, the coupling to planets, Tellus mass).

## 2. Predictions
From the script above, by numerical integration (DOP853) of the equation of motion from the epoch state.

## 3. N's path
N was found near perihelion (2.5 AU from S) and was moving away from the Sun while its sky position moved from lon 52, lat +2 (day 36233) up to lat +46.5 (around day 36375) and then east and south to lon 114, lat +22.5 at the end of the record. Its orbit is steeply inclined (54 deg) and eccentric (0.745).
- Days 36530-36966 (year 100-101.6): it recedes from the Sun (r 3.0 AU at the end of the record) and heads for Jupiter-like P4's orbit (r about 6.2 AU). Its sky position keeps moving east and south (lon 119 -> ~190, lat +21 -> ~0 by day 37100).
- **Day ~36966 (year ~101.2): a very close approach to P4 (0.036 AU, with m_P4 = 1e-3 m_S).** This is the dominant event: the orbit is bent strongly (heliocentric distance grows from ~6 AU to ~10 AU by day 38000) and N is thrown to the outer system; afterwards it moves slowly (r about 10-11 AU through year 120). In the sky it drifts slowly south and east: lon about 214, lat -17 (day 37990), lon 236, lat -31 (39450), lon 275, lat -45 (40910; the southernmost point of the path is lat about -45...-50), then it swings back round to the north (lon 19, lat -21 at 41640; lon 187, lat +2 at 42370) and ends at lon 222, lat -22 (day 43830). The loop in the second half is mostly the parallax of Tellus's orbit against a slowly moving distant body, plus the bend by the encounter. I expect no other close approaches (Tellus >= 1.4 AU, P2 >= 2.1 AU, P3 >= 2.7 AU, P5 >= 10 AU).
- Its motion therefore changes from a fast (up to 0.9 deg/day early on, ~0.4 deg/day at the end of the record) perihelion-like motion to a slow one (a few tenths of a degree per day, retrograde loops each year because of Tellus's own motion).

## 4. Self-assessment
Expected angular errors (my own estimates from a Monte-Carlo over the fit covariance plus the model variants below):
- after 1 year (day ~36890, before the P4 encounter): about 0.1 deg (MC median 0.03 deg, worst case 0.2);
- after 5 years (after the encounter): about 2 deg (MC median 1.0, 84th percentile 2.2, worst 10);
- after 20 years: about 10 deg, and could be several tens of degrees (MC median 5 deg, 84th percentile 11, one branch of the MC gives errors of 100+ deg at year 15 because the outcome of the encounter is chaotic-ish).
The post-encounter path depends on the exact P4 mass, the exact position of P4 (stage-1 ephemeris), and on whether N's coupling to the planets is scaled by (1-beta): I tested the alternative (planet terms multiplied by 1.32); the encounter distance stays 0.035 AU and the predicted path changes by 3 deg at year 5-6 and by about 10 deg at year 20. I took the stage-2 version (not scaled), which has less support from evidence than the rest (I could not test it reliably on other bodies).

What would prove me wrong: (1) N not arriving where predicted at the first year (>0.5 deg off at day 36890 would mean beta or the state is wrong); (2) N's path being deflected at day ~36966 by an amount different from that of a 1.0e-3 m_S body at 0.036 AU (then P4's mass for faint bodies is different, i.e. the planet coupling of faint bodies differs from the Sun coupling); (3) a beta that is not constant (N behaves differently at larger distance, e.g. beyond 6 AU).
