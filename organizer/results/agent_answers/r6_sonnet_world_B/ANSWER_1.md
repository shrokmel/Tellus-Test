# ANSWER_1 - Laws of the Tellus sky

All numbers below are what `predict.py` uses. The model reproduces all 156,000 sightings of all 7 bodies over
100 years at the measurement noise (rms 1.00 arcmin in each angle, no trends, no body worse than 1.01').

## 1. Description of the motions

Frame. The coordinates are tilted relative to the plane in which the Sun moves. The Sun's path is an exact great
circle whose pole is at lon = 90 deg, tilted EPS = 21.7015 +- 0.0007 deg from the frame's pole (this is why the
Sun's latitude swings between -21.70 and +21.70 deg). I rotate about the x axis by EPS to an "ecliptic" frame (call
its coordinates lon_e, lat_e). The plane does not move: its pole drifts by less than 0.01 deg in 100 years.

- **S (Sun)**: moves eastward along the ecliptic, 0.985599 deg/day, period 365.256 d (sidereal year). Speed is
  non-uniform (fastest when at lon_e ~ 273 deg, slowest half a year later): an ellipse of eccentricity 0.0207 seen from
  Tellus. Its ecliptic latitude is 0 within noise; nothing changes over 100 years. (It is really Tellus orbiting S:
  the Sun's apparent direction is the reverse of Tellus's heliocentric direction.)
- **P1**: inner wanderer, never farther than 29 deg from S, swings east and west of S, with retrograde loops near
  inferior conjunction. Orbital period 110.26 d (synodic period 360/(3.265-0.9856) = 158 d).
  Orbit essentially in the ecliptic plane, e = 0.060, a = 0.450.
- **P2**: outer wanderer, period 809.6 d (2.2 yr), retrograde loop at each opposition, e = 0.050, a = 1.70,
  lat_e within +-0.1 deg (orbit in the ecliptic plane, inclination consistent with 0).
- **P3, P4, P5**: outer wanderers with periods 2190 d (6.0 yr), 5640 d (15.4 yr), 14240 d (39 yr), each with a
  retrograde loop at opposition (loop shorter in time and bigger the nearer the planet). Orbit inclinations 1.50, 1.00,
  2.00 deg to the ecliptic, so their latitudes reach +-2.2, +-1.2, +-2.2 deg. Eccentricities 0.040, 0.050, 0.040.
- **M** (large bright body): goes round the sky in 27.3 d (13.19 deg/day) with lat_e up to +-5.2 deg, i.e. an orbit
  around Tellus (not around S) tilted 5.00 deg, e = 0.055. Its orbit is strongly perturbed by S: the line of nodes
  regresses (moves westward) by 0.05295 deg/day = one turn per 18.61 yr, the perihelion advances too. These
  show as slow drifts over the 100 years.
- Changes over the 100 years: no secular change of the Sun's orbit or of the frame. The planets show slow
  perturbations of their elliptical orbits (apsides move by fractions of a degree to ~2 deg per century, mean
  motions wobble by a few arcmin) caused by their mutual attraction (largest for P3, P4, P5 which interact with each
  other, P4 being the heaviest). M's node/apse motion is the large one.

## 2. Laws as equations

**Law of gravity (one law for everything).** Every body i accelerates toward every other body j:

    d^2 r_i / dt^2 = sum_{j != i} GM_j (r_j - r_i) / |r_j - r_i|^3 ,    GM_j = GM_S m_j

with r in AU, t in days, GM_S = 2.9592e-4 AU^3/day^2 (equivalently n^2 a^3 = GM_S (1+m): Kepler's third law;
checked: 2-body period T = 2 pi sqrt(a^3/GM_S) matches all planets to <1e-4). The exponent of r (-2) was not fitted
separately; it was assumed, and with it every body fits at noise level over 100 years (a different exponent would
give apsidal precession of every orbit, not seen). No other force (no drag, no extra precession term) is needed.
GM_S is measured only through the unit of length: a_P = 1 AU is fixed by GM_S = 2.9592e-4 (+-1e-6 relative; the
fitted semi-major axes come out as round numbers 0.45, 1.70, 3.30, 6.20, 11.50 in these units).

**Bodies.** S = massive centre (m = 1). Tellus E (m = 3.0e-6, assumed; not constrained by directions alone), M orbits
E (m_M = 3.7e-8, assumed). P1..P5 planets orbiting S. The motion of S about the barycentre is included.

**Observed direction** of body X from Tellus: unit vector of (r_X - r_E), then lon_e = atan2(y,x), lat_e = asin(z/r)
(no light-time, no aberration needed at the 1' level), and the frame transformation to the catalogue frame

    x = x_e ,  y = cos(EPS) y_e + sin(EPS) z_e ,  z = -sin(EPS) y_e + cos(EPS) z_e ,   EPS = 21.7015 deg.

**Initial conditions at day 0** (osculating Kepler elements, heliocentric for E, P1..P5; geocentric for M; a in AU, angles
in deg; varpi = longitude of perihelion, L0 = mean longitude at day 0, i, Omega inclination and node to the ecliptic):

| body | a | e | varpi | L0 | i | Omega | mass m / M_S |
|---|---|---|---|---|---|---|---|
| E (Tellus) | 0.999319 +- 1e-5 | 0.020686 +- 5e-7 | 93.462 +- 0.002 | 246.1448 +- 0.0002 | 0.005 (unconstrained, <0.01) | 30 (undefined) | 3.0e-6 (fixed) |
| M (about E) | 0.0025700 +- 1e-8 | 0.055001 +- 1e-6 | 85.284 +- 0.002 | 283.1385 +- 0.002 | 5.0010 +- 0.0002 | 298.013 +- 0.003 | 3.7e-8 (fixed) |
| P1 | 0.449999 +- 5e-6 | 0.060010 +- 1e-6 | 166.419 +- 0.003 | 218.3185 +- 0.001 | 0 (<0.01) | - | <1e-5 (used 2e-6) |
| P2 | 1.699994 +- 3e-5 | 0.050003 +- 1e-6 | 281.186 +- 0.002 | 114.9178 +- 0.0003 | 0 (<0.01) | - | <1e-5 (used 3e-7) |
| P3 | 3.300001 +- 2e-6 | 0.039993 +- 1e-6 | 42.543 +- 0.001 | 151.6885 +- 0.0003 | 1.4990 +- 0.0001 | 118.667 +- 0.003 | (1.4956 +- 0.0004)e-4 |
| P4 | 6.19994 +- 1e-5 | 0.050009 +- 1e-6 | 211.514 +- 0.002 | 356.6341 +- 0.0003 | 0.9985 +- 0.0001 | 163.244 +- 0.003 | (1.0000 +- 0.0001)e-3 |
| P5 | 11.4998 +- 1e-4 | 0.040006 +- 1e-6 | 167.698 +- 0.003 | 77.8088 +- 0.0003 | 2.0005 +- 0.0001 | 73.202 +- 0.002 | (2.498 +- 0.001)e-4 |

(Uncertainties are formal 1-sigma from the least-squares covariance; real uncertainties are probably 2-3x larger
because of degeneracies. They are the osculating elements at day 0 of the N-body integration, not mean elements.)
Time is in days from day 0; "day" t is exact integer days; the integrator is a 6th-order symplectic Yoshida scheme with
8 steps/day (in tests, 4 steps/day differs from 8 by 0.09 arcmin for M and <0.001' for the others).

Not independent laws but consequences: Kepler's laws for each body, M's nodal regression of 18.6 yr, the slow
apsidal/mean-motion changes of the planets (planet-planet attraction; the mass ratios show the heavy planet is P4,
~1/1000 of S).

## 3. Predictions

`predict.py` (run with the given python) integrates the above equations from the day-0 state to day 43815 and writes
`predictions.csv` (lon, lat in the catalogue frame, degrees). It takes ~1 min. No observation file is read; the
constants are the ones in section 2.

## 4. Self-assessment

- Inverse-square Newtonian gravity applies to all bodies, including M around Tellus and the planets among
  themselves: 97%. (Everything fits to noise over 100 years with 1 law; a hidden extra massive body beyond P5 or very
  small corrections to the exponent would be the thing I might have missed.)
- Sun/Tellus orbit (ellipse e = 0.0207, plane fixed, EPS = 21.7015): 99%.
- Planet orbits and elements P1, P2 (nearly circular, in-plane): 98%; P3-P5 elements and masses: 95%.
- P4 mass 1.000e-3 and P5/P3 masses: 90% (well measured via perturbations; P1, P2, Tellus, M masses are NOT measured
  by these data: only upper limits ~1e-5).
- Predictions for years 101-120: expect errors of order the noise (1 arcmin) for S, M, P1-P4; P3, P4, P5 maybe 1-3'
  because small element errors grow over 20 years. M is the one most sensitive to the numerical scheme.

What would prove me wrong: a sighting of P5, P4 or P3 that deviates from the prediction by several arcminutes
growing with time (would mean a missing body or a change of the force law); a significant latitude for P1 or P2
(inclination I set to 0); Sun's latitude departing from 0 (Tellus orbit tilting); a change of the node/precession
of M inconsistent with 18.6 yr regression; a new body appearing in the data; any drift of EPS larger than 0.01 deg.
