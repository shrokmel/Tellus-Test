# ANSWER 1

Frame: directions are ecliptic longitude/latitude in the data frame (pole N = [-3.138e-5, 0.369770, 0.929123], normalised). Time t is days since observations began. Angles in degrees unless stated.

## 1. Description of the motions

- **S (Sun, the body that makes daylight).** Its apparent path is a great circle, inclined 21.70 deg to the frame's pole, with the node near lon 90 deg. It is a clean Kepler orbit: period 365.2576 d, eccentricity 0.01998. Residuals over all 100 years are 61" rms, which is the 1' measurement noise. No change over 100 years is detectable.
- **P1 (inner, planar).** Coplanar with the S-plane (geocentric latitude scatter is noise only). Kepler, period 110.2578 d, eccentricity about 0.0601. Over 100 years the period is constant to about 3e-4 d (window-to-window scatter), the perihelion direction does not measurably move, and the eccentricity drifts slowly and steadily downward, from 0.06009 to 0.05991. Residuals are at noise level (about 60" rms) in every 1500-day window.
- **P2 (outer, nearly planar).** Geocentric latitude is noise-level, so the orbit is essentially in the S-plane. Period about 809.7 d, eccentricity about 0.050. Fits to about 135-150" rms, with a systematic offset of about +110" in most decades. The model is not complete.
- **P3 (outer, inclined about 1.5 deg).** Period about 2189 d, eccentricity about 0.042. A static Kepler model misfits it badly (about 690" rms, with decade-mean offsets of several hundred arcsec). Its residuals show structure that is not captured by a fixed orbit.
- **P4 (outer, inclined about 1.0 deg).** Period about 5640 d (about 15.4 yr), eccentricity about 0.049. Fits to about 165" rms; residuals vary by decade by about +/-50".
- **P5 (outer, inclined about 2.0 deg).** Period about 14169 d (about 38.8 yr), eccentricity about 0.041. Fits to about 265" rms, with large decade-to-decade offsets (up to about 450").
- **M (large bright body, orbits Tellus).** Inclination about 5 deg to the S-plane, sidereal period about 27.29 d (fitted on the full dataset as 27.2897 d; the latitude periodogram gives about 27.18-27.29 d). Its node regresses by roughly 20 deg per year and its perigee advances on a similar timescale. Fixed-element fits fail over any window longer than about 60 days. A single static Kepler orbit does not describe M over the 100 years, and 60-day windows still leave about 0.6-0.8 deg rms with some failed fits. I could not model M's secular motion, so its element changes are unresolved.

Over 100 years: S and P1 are steady apart from P1's slow eccentricity decrease. P2-P5 show slowly changing mean offsets that I did not model. M's orbit geometry changes continuously and is the largest unresolved effect.

## 2. Laws (equations and constants)

Notation: t in days; M(t) = mean anomaly; E = eccentric anomaly from Kepler's equation M = E - e sin E; nu = true anomaly; r = a(1 - e cos E); body directions are normalised vectors. Frame basis: EX = normalised (z - (z.N)N), EY = N x EX.

**Law 1 (Kepler's laws for Tellus / Sun).** Tellus moves on a Kepler orbit about S; the geocentric Sun has longitude = varpi_S + nu_S(t) and latitude determined by the S-plane.
- P_E = 365.257550 d (uncertainty not formally estimated)
- e_E = 0.019980
- varpi_S = 3.9069 deg
- M0_S = 152.2598 deg
- a_E = 1 AU (unit length). Residual rms 61" (noise).

**Law 2 (Kepler's third law, Newtonian gravity about S).** For each planet a_p = (P_p/P_E)^(2/3) AU, so the geometry is fixed by its period.

**Law 3 (planet P1, planar orbit).** Heliocentric position r = a_1 (cos u, sin u, 0) in the S-plane, u = varpi_1 + nu_1, with P_1 = 110.2578 d, e_1 = 0.06009 (decreasing; slope about -1.8e-6 per year), varpi_1 = 256.468 deg, M0_1 = 51.841 deg. Geocentric direction = (r_1 - r_Tellus) normalised. Residual 58" rms.

**Laws 4-7 (planets P2-P5, 3D Kepler about S).** Heliocentric Kepler orbits with (P, e, i, Omega, omega, M0):
- P2: P = 809.67 d, e = 0.0501, i = 0.006 deg, Omega = 306.279 deg, omega = 65.598 deg, M0 = 193.013 deg. Residual about 150" rms.
- P3: P = 2189.15 d, e = 0.0421, i = 1.493 deg, Omega = 208.253 deg, omega = 284.133 deg, M0 = 109.247 deg. Residual about 690" rms (poor).
- P4: P = 5639.82 d, e = 0.0493, i = 1.003 deg, Omega = 253.056 deg, omega = 49.059 deg, M0 = 144.474 deg. Residual about 165" rms.
- P5: P = 14169.46 d, e = 0.0413, i = 1.999 deg, Omega = 162.952 deg, omega = 98.863 deg, M0 = 266.175 deg. Residual about 265" rms.
Uncertainties on these fitted constants are not formally estimated; from the decade-to-decade scatter they are of order 1e-3 in e and 0.1 deg in angles, and the period errors are larger for P3-P5.

**Law 8 (Moon M, Kepler orbit about Tellus).** Geocentric orbit with P_M = 27.3302 d, e = 0.0647, i = 4.927 deg, Omega = 281.24 deg, omega = 291.33 deg, M0 = 147.04 deg. These are the elements of the last 60-day window and are low confidence: the node and perigee are drifting, and the elements do not extrapolate.

Gravitational law: all bodies obey Newton's law of gravitation about S (inverse square); the planets' mutual perturbations and the Sun's perturbation of M are the likely source of the unmodelled residuals. I did not fit masses or an N-body integration.

## 3. Predictions

`predictions.csv` (days 36555 onward, 1701 rows) was generated by `predict.py` in this folder, which implements Laws 1-8 with the constants above.

Check against the last observed days (days 36000-36400, in-sample, so not an out-of-sample test), median angular error:
- S 0.025 deg; P1 0.022 deg; P4 0.020 deg; P5 0.030 deg; P2 0.086 deg; P3 0.24 deg; M 2.7 deg.
Predictions for M will degrade quickly beyond the sampled window (the Moon's node and perigee move about 20 deg per year). Predictions for P3 and P2 are less reliable than their short-window accuracy suggests, because their residuals show secular structure.

## 4. Confidence and falsifiers

- Law 1 (Sun, Kepler, P_E = 365.2576 d): 95%. Would be falsified by a residual trend beyond about 60" over decades, or a change in the S-plane tilt.
- Law 3 (P1 planar Kepler, P = 110.2578 d): 90%. Falsified by out-of-plane geocentric latitude for P1 (none seen above noise), or a period change. The eccentricity decrease is about 2 sigma in per-window scatter and may be real or an artefact; 50% that it is a real secular trend.
- Laws 4-7 (P2, P4 3D Kepler): 60%. Residuals are 2-3 times the noise and have systematic decade structure; the static orbits are approximations. P5: 40%. P3: 25%.
- Law 8 (Moon, Kepler about Tellus with fixed elements): 15%. The fixed orbit is known to be wrong over long spans; the true motion is perturbed by S and changes on the time scale of a few months.
- Overall gravitation (inverse-square, Newtonian N-body): untested here; I did not fit an N-body model, so this is a hypothesis, not a result.

What would prove me wrong: a P1 geocentric latitude out of the S-plane beyond 0.1 deg; a new decade-scale residual in the Sun beyond 60" rms; or out-of-sample predictions (year 101-120) for P2-P5 and M that miss by more than the stated residuals.
