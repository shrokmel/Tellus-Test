# Answer: the sky of Tellus

All numbers come from least-squares fits to the 156,000 observations (100 years). The final model reproduces every body, in both lon and lat, to the 1 arcminute measurement noise (rms 0.99-1.02 arcmin per coordinate for each of S, P1-P5, M over the full 36,525 days, with no trend in the last decade). `predict.py` implements exactly the equations below; `work/` holds the analysis scripts.

## 1. Description of the motions

* **Reference plane.** S (the Sun-like body) always stays in one great circle tilted by **eps = 21.702 deg** to the equator of the (lon, lat) frame, with the node at lon = 0/180. Its ecliptic latitude is zero within noise (rms 1.0 arcmin) and the tilt does not change over 100 years. P1 and P2 also lie in this plane (inclinations 0 within about 1e-5 rad); P3, P4, P5 are inclined by 1.5, 1.0 and 2.0 deg, with fixed nodes. M reaches |lat| = 26.85 deg, i.e. its path is inclined about 5.0 deg to the ecliptic.
* **S**: moves eastward, 0.9856 deg/day (period 365.26 d), unevenly (fast in one season, slow in the opposite): an ellipse seen from Tellus. This is Tellus orbiting S, reflected: S's apparent lon is the opposite of Tellus's heliocentric longitude. Eccentricity 0.0200. Its apsidal line barely moves.
* **P1-P5**: wanderers that go around S, not Tellus. Seen from Tellus they show direct motion, stations and retrograde loops. P1 (period 0.302 yr = 110 d) is an inner body that oscillates around S (never far away in angle, like an inner planet), P2 (2.22 yr) is an outer body with retrograde loops once per synodic period, P3 (5.99 yr), P4 (15.4 yr), P5 (39.0 yr). Retrograde loops happen when Tellus overtakes them, which is the sign of heliocentric orbits.
* **Slow changes over 100 years**: a pure Kepler ellipse per body fits the first ~4 years to noise, but over the century P3 drifts by up to 0.5 deg, P5 and P4 by tens of arcmin, P2 by a few arcmin. The residuals have periods equal to the synodic combinations of the planets (half the P3-P4 synodic period 4.9 yr, P3-P5 7.2 yr, P4-P5 26 yr). So the planets perturb each other gravitationally (P4 is by far the strongest perturber). There is no sign of any other secular change: no change of the obliquity, no change of orbital planes beyond what the mutual attraction produces, and no extra precession.
* **M**: a large bright body going around the sky fast, 13.18 deg/day (27.3 d), always near S's plane, with a large irregular motion (the apparent speed varies; its path precesses). It is a satellite of Tellus (its geocentric orbit is perturbed by S: the line of nodes regresses and the apsides advance, an 'evection'-like irregular motion). It is not seen to affect S or the others.
* **Visibility** gaps (within 15 deg of S, clouds) are just sampling.

## 2. Laws (equations and constants)

**Law 0 (frame).** Define the ecliptic frame from the observed frame by a rotation about the x axis (lon = 0 direction) by eps:

    (x_e, y_e, z_e) = (x, y cos eps - z sin eps, y sin eps + z cos eps),   eps = 21.702 +- 0.001 deg
    (x, y, z) = (cos lat cos lon, cos lat sin lon, sin lat)

The fit of S's latitude in this frame is zero (|mean| < 0.01 arcmin, rms 1.0 arcmin, constant for 100 years).

**Law 1 (gravity).** Every body i attracts every body j with Newtonian gravity:

    d^2 r_j / dt^2 = sum_{i != j} G m_i (r_i - r_j) / |r_i - r_j|^3

with point masses, in an inertial frame fixed against the stars. Only directions are observed, so lengths are relative: the unit of length is Tellus's orbital semi-major axis, the unit of time is the day, and the central body S is the unit of mass. Then G m_S = n_T^2, with n_T given below (this is Kepler's third law, n^2 a^3 = G m_S (1 + m)).

The direction of body B seen from Tellus is the direction of (r_B - r_T), rotated to (lon, lat) with Law 0.

**Law 2 (initial conditions at day 0).** The bodies start on osculating heliocentric Kepler ellipses given by n (rad/day), mean longitude L0 at day 0, eccentricity e, longitude of periapsis varpi, node Omega and inclination i (to the ecliptic). Semi-major axis a = (G m_S (1+m) / n^2)^(1/3) (so a = (n_T/n)^(2/3) for each planet in units of Tellus's). The constants (1-sigma from the covariance of the fit):

| body | n (rad/day) | period | L0 (rad) | e | varpi (rad) | Omega (rad) | i (deg) | a (Tellus=1) | m / m_S |
|---|---|---|---|---|---|---|---|---|---|
| Tellus (S seen from it) | 0.017202503 +- 1e-9 | 365.261 d | 4.29648 +- 4e-6 | 0.020001 +- 1.3e-6 | 1.63693 +- 1.1e-4 | - | 0 (defines plane) | 1 | not constrained (taken 0) |
| P1 | 0.05698639 +- 4.5e-9 | 110.26 d | 3.81040 +- 3e-5 | 0.059999 +- 5e-6 | 2.90523 +- 8e-5 | - (i = 5e-6 deg-ish, taken as i cos Omega = -8.8e-6, i sin Omega = -8.5e-6 rad) | ~0.0005 | 0.4500 | test particle (0) |
| P2 | 0.007760953 +- 1.1e-9 | 2.2166 yr | 2.00571 +- 4e-6 | 0.049990 +- 1.4e-6 | 4.90802 +- 4e-5 | n/a | 0.0012 +- 0.0001 (~0) | 1.6998 | (-0.5 +- 1.5)e-7 (consistent with 0) |
| P3 | 0.002869830 +- 3e-9 | 5.994 yr | 2.64744 +- 4e-6 | 0.040006 +- 1.7e-6 | 7.02555 +- 3.4e-5 | 2.07084 +- 1e-4 | 1.4991 +- 0.00014 | 3.2995 | (1.5146 +- 0.0005)e-4 |
| P4 | 0.001114855 +- 1.9e-9 | 15.43 yr | -0.05872 +- 4e-6 | 0.049992 +- 2.3e-6 | 3.69154 +- 3.6e-5 | 2.84887 +- 1.5e-4 | 0.99813 +- 0.00015 | 6.198 | (1.0003 +- 0.0005)e-3 |
| P5 | 0.000441147 +- 1.4e-9 | 38.99 yr | 1.35800 +- 4.6e-6 | 0.039997 +- 1.7e-6 | 2.92606 +- 5.8e-5 | 1.27743 +- 7.5e-5 | 2.0005 +- 0.00015 | 11.50 | (2.491 +- 0.007)e-4 |

(Angles in rad; the planets' elements are heliocentric ecliptic, epoch day 0. Some of these hint at round values: e = 0.02, 0.06, 0.05, 0.04, 0.05, 0.04; i = 0, 0, 1.5, 1.0, 2.0 deg; a = 0.45, 1.70, 3.30, 6.20, 11.50; m = 1.5e-4, 1.0e-3, 2.5e-4. These are what the fit gives within errors.)

The fitted P1 inclination from the P1 fit is not zero in a strict sense (about 1e-5 rad), but it is negligible for prediction.

**Law 3 (M is Tellus's satellite).** M moves under Law 1 in the field of Tellus plus S's attraction. In the geocentric frame (r = position of M relative to Tellus, s = position of S relative to Tellus):

    d^2 r/dt^2 = - mu r / |r|^3 + G m_S [ (s - r)/|s - r|^3 - s/|s|^3 ],   mu = G m_T (+ m_M)

The tidal term depends on the ratio (n_T/n_M)^2, so the angles do not depend on the (unobservable) scale; mu = 3.0e-6 G m_S was used as a convention. Osculating elements at day 0 (geocentric, in the ecliptic frame): n = 0.2300980 +- 2.6e-8 rad/day (osculating; the mean sidereal period is 27.3 d), L0 = 4.94174 +- 4e-6, e = 0.054992 +- 1.4e-6, varpi = 7.77167 +- 2.6e-5 rad, Omega = 2.05887 +- 3e-5, i = -0.087259 +- 2.2e-6 rad (= 5.000 deg with sign convention). Nothing beyond Newtonian gravity (no oblateness, no drag) was needed: the residual is exactly the noise.

**Integration.** Laws 1-3 are integrated numerically (DOP853, rtol 1e-10) in `predict.py`: first the six-body system S, Tellus, P2-P5 (barycentric), then P1 as a test particle in that field, then M. The solution is the 'equation'; no closed form exists for the perturbations, but the unperturbed limit of each orbit is the Kepler ellipse r = a(1 - e cos E), M = E - e sin E, M = n t + L0 - varpi.

Fit statistics: all 7 bodies at 0.99 - 1.02 arcmin rms. The S, P2-P5 system was fit jointly (32 parameters), P1 and M separately given the joint ephemeris.

## 3. Predictions

`predictions.csv` (days 36555-43815, every 30 days, all 7 bodies) is produced by `python predict.py` from the constants above.

## 4. Self-assessment

| Claim | Confidence |
|---|---|
| Bodies P1-P5 orbit S (not Tellus) on Keplerian ellipses; Tellus orbits S | 99% |
| All orbits are coplanar to within 2 deg, the reference plane is tilted 21.70 deg to the frame and fixed | 98% |
| Inverse-square Newtonian gravity incl. planet-planet attraction (Kepler's 3rd law holds for all) | 95% (the P3-P4-P5 perturbation pattern and masses are fitted with errors of a few 1e-3 relative; an exponent other than 2 or a slightly different law was not tried, the test is the zero-parameter match of n^2 a^3 for P1-P4 and the joint fit at the noise level) |
| Masses: P4 ~ 1.0e-3, P5 ~ 2.5e-4, P3 ~ 1.5e-4, P2 ~ 0 | 90% (P2 could be tiny but non-zero; P1's and Tellus's masses are not measurable) |
| M is a satellite of Tellus perturbed by S (Newtonian) | 90% |
| Extrapolation to years 101-120 is accurate to a few arcmin | 85% (the errors grow slowly from the parameter errors, estimated < 1 arcmin for most bodies, a few arcmin for P3-P5 because of the mass and phase errors) |
| Small effects not modelled: P1's mass, Tellus's mass and the Tellus-M barycentre wobble of S (S fit is perfect at 1 arcmin so these are below ~0.3 arcmin) | - |

What would prove me wrong:
* A systematic drift of S's ecliptic latitude, or of the orbital planes, beyond mutual perturbations (e.g. changing obliquity) in years 101-120.
* A prediction error for P3, P4 or P5 growing faster than the parameter covariance implies (a missing body or non-inverse-square force), or a new periodic signal in P3's residuals at a synodic period not generated by P4/P5.
* M departing from the tidal model: for instance a regular difference in the node regression rate would mean extra forces (oblate Tellus) or a different law.
* P2's mass: if later data show a P2-induced perturbation on P3/P4 the zero mass is wrong.
