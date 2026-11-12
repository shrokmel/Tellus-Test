# ANSWER_1

## 1. Description of the motions

- The frame is equatorial (fixed to the stars). The Sun S moves on a near-circular ellipse (e = 0.020) with period 365.26 d. Its orbital plane (the "ecliptic") is tilted 21.7 deg to the frame, so S has latitude of +/-21.7 deg. The Sun's orbit does not drift over the 100 years.
- M (Moon) goes round Tellus every 27.29 d (sidereal). Its geocentric path is strongly perturbed by the Sun (lon residual from pure Kepler about 3.4 deg rms). Its node regresses (latitude extremes move between about 19.6 and 26.9 deg).
- Inner body P1 never gets more than 29 deg from S (max elongation). Heliocentric period 110.26 d, a = 0.450 AU, e = 0.060.
- P2: period 809.7 d, a = 1.700 AU, e = 0.050.
- P3: period 2189 d, a = 3.299 AU, e = 0.042. It is the only planet with a systematic residual. Its residual has its main harmonic at 2(n3 - n4), where n3 and n4 are the mean motions of P3 and P4. That points to a gravitational perturbation by P4.
- P4: period 5640 d (15.4 yr), a = 6.203 AU, e = 0.049.
- P5: period 14169 d (38.8 yr), a = 11.466 AU, e = 0.041.
- All of these orbits are nearly coplanar with the Earth-Sun plane (inclinations below 0.4 rad with respect to the equator frame, i.e. about the 21.7 deg tilt of the ecliptic).
- Over the 100 years the elements are constant within the noise. Window fits (5000 to 10000 days) show: P1 and P5 are constant to the noise (about 0.8 arcmin rms). P2 shows a small drift of omega (+0.02 rad over 100 yr) and of the node (+0.0004 rad), which may be real. P3's e and omega wobble by about +/-0.001 and +/-0.04 rad, consistent with the P4 perturbation.

## 2. Laws (equations and constants)

Units: angles in radians, time in days, distance in AU. Positions are heliocentric unless stated.

**Law 1 (gravitation with one common GM for the Sun).** The Sun's gravitational parameter fixes every orbital period:
  GM_S = n_E^2 * (1 AU)^3,  with n_E = 0.0172021 rad/d (so GM_S = 2.9598e-4 AU^3/d^2).
  Kepler's third law holds for all bodies with this GM: a^3 n^2 / GM_S = 1.000 +/- 0.002 for P1, P2, P3, P4, P5 (checked from the fits).

**Law 2 (Kepler orbits about the Sun).** For each body the heliocentric position is
  r(t) = Kepler(a, e, i, Omega, omega, M0, n; t),  with M = M0 + n t, and standard Kepler/ perifocal rotation.

Constants (fitted; uncertainties are the spread between 5000 or 10000 day windows, or the Earth-fit residual):

| body | a [AU] | e | i [rad] | Omega [rad] | omega [rad] | M0 [rad] at day 0 | n [rad/d] | T [d] |
|---|---|---|---|---|---|---|---|---|
| Tellus (Earth) | 1 (fixed) | 0.01998 +/- 0.0002 | 0.37877 +/- 0.0001 | 3.14169 +/- 0.0002 | -1.50264 +/- 0.001 | 15.22383 | 0.0172021 +/- 0.0000003 | 365.258 |
| P1 | 0.45000 +/- 0.00002 | 0.05999 +/- 0.0001 | 0.37876 +/- 0.00005 | 3.14163 +/- 0.0001 | -0.23580 +/- 0.0003 | 7.18774 +/- 0.0002 | 0.0569863 +/- 0.0000005 | 110.258 |
| P2 | 1.70006 +/- 0.0003 | 0.05005 +/- 0.0001 | 0.37885 +/- 0.0001 | 3.14184 +/- 0.0002 | 1.77783 +/- 0.001 | 3.36880 +/- 0.0003 | 0.0077602 +/- 0.000001 | 809.67 |
| P3 | 3.29941 +/- 0.0003 | 0.04206 +/- 0.001 | 0.39175 +/- 0.0002 | 3.08151 +/- 0.0003 | 3.93707 +/- 0.03 | 1.90677 +/- 0.03 | 0.0028701 +/- 0.000001 | 2189.1 |
| P4 | 6.20316 +/- 0.0007 | 0.04926 +/- 0.0002 | 0.39555 +/- 0.0001 | 3.12841 +/- 0.0002 | 0.57280 +/- 0.002 | 2.52157 +/- 0.002 | 0.0011141 +/- 0.0000003 | 5640 |
| P5 | 11.46612 +/- 0.001 | 0.04125 +/- 0.0002 | 0.36997 +/- 0.0001 | 3.04934 +/- 0.0001 | 6.22644 +/- 0.003 | 4.64547 +/- 0.003 | 0.0004434 +/- 0.0000003 | 14169 |

Law 2 is the Kepler form of Laws 1 and 2 combined; the Earth row is the Sun's geocentric orbit.

**Law 3 (observation geometry).** The Sun's geocentric direction is -r_E / |r_E|. A planet's geocentric direction is (r_p - r_E)/|r_p - r_E|. Only directions are observed.

**Law 4 (Moon).** The Moon's geocentric direction is r_M / |r_M| where r_M follows a geocentric Kepler orbit (e = 0.00100, i = 0.38259, Omega = 3.15129, omega = 5.60697, M0 = 2.48693, n = 0.2302405 rad/d, T = 27.2897 d, in a scale-free form), plus an empirical harmonic correction:
  delta_lon * cos(lat) + i delta_lat = sum_k [ c_k^+ exp(+i w_k t) + c_k^- exp(-i w_k t) ],
with 39 frequencies w_k built from the Moon's mean motion n_M, the Sun's n_S = n_E, the lunar node (T = 6798 d) and perigee (T = 3233 d). The coefficients c_k are in `predict.py`. The harmonic terms are the Sun's perturbation on the Moon (evection, variation, annual equation). The in-sample residual is 18 arcmin rms (per component about 14 arcmin lon, 21 arcmin lat).

**Law 5 (P3 perturbation).** P3's Kepler orbit is corrected by a harmonic set at frequencies n3 - n4, 2(n3 - n4), etc. (9 frequencies, 18 coefficients in `predict.py`). This is an empirical stand-in for the P4 perturbation.

I did not establish a full N-body model with fitted masses (see section 4). The planet masses are not determined.

## 3. Predictions

`predictions.csv` (243 days, years 101-120, all seven bodies) is produced by `predict.py` using only Laws 1-5 with the constants above and in the script. Predicted lon/lat in degrees in the equatorial frame.

In-sample check against the observations (rms, arcmin, lon / lat): S 1.0 / 1.0; P1 1.1 / 1.0; P2 2.9 / 1.2; P3 8.2 / 2.6; P4 3.1 / 1.3; P5 4.8 / 1.6; M 14 / 21.

Expected prediction error (extrapolation over 1 to 20 years beyond the data): about 1-3 arcmin for S, P1, P4 and P5; about 10 arcmin for P2 and P3; about 0.5 to 1 degree for M (the harmonic set was not validated out of sample).

## 4. Confidence and what would prove me wrong

- **Law 1 / 2 (Kepler orbits, common Sun GM), about 90%.** Each planet fits at the noise level in its own 5000-10000 day windows, and a^3 n^2 / GM_S = 1 within 0.2% for all bodies. Proof of error: a planet's a^3/T^2 ratio differs from the Sun's by more than 0.5%, or a 1-arcmin-level residual appears in a planet with no perturber.
- **Constancy of the orbits over 100 years, about 75%.** Elements are constant within noise for S, P1, P4, P5, and roughly for P2 (small omega/node drift). Proof of error: an omega or node drift for P1 or P5 larger than 0.001 rad per century, or a systematic change in the Sun's period.
- **P3 perturbation by P4 (harmonic model), about 60%.** The residual frequency matches 2(n3 - n4) and its amplitude is about 10 arcmin, but I did not fit a full N-body model, so the perturbing mass is not determined. Proof of error: P3 residuals that do not shrink when the 2(n3 - n4) harmonic is included, or a predicted P3 direction that misses by more than 10 arcmin.
- **Moon harmonic model, about 50%.** The fit uses 39 frequencies, including node and perigee periods I took from the Moon's known cycles (18.6 yr and 8.85 yr) rather than fitting them. The model reproduces the in-sample data to 18 arcmin, but I did not check it out of sample. Proof of error: Moon residuals larger than about 30 arcmin in the predicted years, or a drift in the node or perigee periods.
- **Overall**: the most likely failure is the Moon prediction. The planets' Kepler laws are the most reliable.

Physically consistent N-body integration with fitted masses was not completed: the restricted Moon integration was too slow (about 15 s per evaluation) and the least-squares fit did not move from its start. The harmonic corrections above are empirical substitutes for that.
