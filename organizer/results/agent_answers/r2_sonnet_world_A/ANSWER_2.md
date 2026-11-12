# Stage 2: the faint bodies

Files: `predict_2.py` (run it from this folder; it imports the stage-1 script `predict.py`) writes `predictions_2.csv` (19 bodies x 243 days). Analysis scripts are in `work/`.
Honest summary up front: the stage-1 laws are untouched and still explain S, P1-P5 and M to the noise. The faint bodies C1-C12 are test particles moving under the same gravity, but (i) each one needs its own factor on the attraction of S, (ii) long-term orbits change at planetary encounters, and (iii) one close Earth flyby (C1, and probably C7) is NOT explained by the Earth mass that the other bodies prefer. Some bodies (C2, C4, C5, C8, C9, C12) are fitted to noise level on their recent history; C1, C3, C6, C7, C10 and C11 are only partly constrained and their predictions are much less reliable.

## 1. Test of the stage-1 laws

**Holds (stage-1 bodies).** The stage-1 model (six-body Newtonian system S, Tellus, P2-P5 + P1 and M as test particles) reproduces the 100-year record of S, P1-P5, M with rms 0.99-1.02 arcmin per coordinate, unchanged. Adding a Tellus mass of 3e-6 m_S to the stage-1 N-body and refitting the elements leaves the rms at 0.9995, so stage 1 cannot see Tellus's mass (the elements absorb it).

**Holds for the faint bodies (with additions).** A faint body is a massless test particle that obeys inverse-square gravity from S plus the planets. With the stage-1 planets (P3, P4, P5 masses as in stage 1) and a Tellus mass of 3e-6 m_S, a single heliocentric orbit fits a long run of observations at the noise level (median vector residual ~1.1-1.5 arcmin; noise alone gives 1.17):

| body | data fitted to noise level | median residual (arcmin) | comment |
|---|---|---|---|
| C2 | all 10 apparitions, days 6456-36525 (3 apparitions are 2-24 arcmin) | 1.6-2.4 (two apparitions 16 and 24) | a=4.496, e=0.446, i=21.64 deg, beta ~ 0 |
| C4 | last ~9000 days (35 chunks in total) | 1.2-1.5 | near-Earth orbit, degrades to 2-5 arcmin over 25,000 d and to degrees earlier |
| C5 | days 20000-34300 (5 apparitions) | 1.3-1.6 | the orbit before ~day 20000 is different (planetary encounter) |
| C8 | days 35872-36402 (2 apparitions) | 1.9-2.9 | |
| C9 | days 20000-34200 (4 apparitions) | 3-6 | slightly above noise |
| C12 | days 25143-34857 (2 apparitions) | 1.1-1.3 | earlier apparitions fail (hundreds of arcmin) |
| C10 | 3 apparitions | 6.8, 12.7, 13.9 | not at noise |
| C11 | last apparition only (35516-35741) | 1.2-2.0 | 6 earlier apparitions fit at 2-30 arcmin, last 2 apparitions are mutually inconsistent at 10-15 arcmin |
| C3 | 3 apparitions (days 5000-14300) | 5, 11, 31 | no data after day 14297 |
| C6 | last apparition only (34575-34907) | 1.2 | the 2 earlier apparitions are not reproduced |
| C7 | last chunk only (~day 35800-36420) | 5.0 | earlier chunks fit at 1.2-1.5 but with a different orbit (a changes by 3 % at a 0.068 AU Earth flyby on day ~35777) |
| C1 | only after the flyby (days 8200-8331) | 1.3 | see below |

**Fails / forced changes.**
1. *Pure Newton from S fails for some bodies.* C11 (a 520-day apparition) cannot be fitted with 1/r^2 attraction of S at full strength: rms 44 arcmin per point; with a radial extra acceleration it fits at 1.4 (radial A = -46.8e-6 AU^3/d^2 / r^2 = 0.158 GM_S/r^2, tangential and normal components 0.03 and 0.007, i.e. zero). Fitting the exponent gives r^-1.9998 +- small, so the extra term is also inverse-square: the effective attraction of S is (1-beta) GM_S with beta_C11 = -0.150. Independent apparition pairs give beta_C11 = -0.149...-0.155. Other bodies have beta = 0 within 1e-3 (C2, C4, C5, C7), or clearly not 0 (C3 -0.25, C8 +0.297, C9 +0.082, C12 -0.296). This is the first law that stage 1 did not need.
2. *A body's orbit is not constant over the century.* C2, C5, C11, C12, C7 and C8 require different orbits before and after certain dates (e.g. C5 before day ~20000; C12 before 25000; C8 before 27000); 2-body fits of single apparition pairs are often at noise while the N-body orbit continued through an encounter does not match. I attribute this to planetary encounters (P3/P4 for C5/C11, Tellus for C7, C1), but I could not prove it for each case.
3. *Earth flyby of C1 (and C7) is not explained.* C1's direction jumps by 144 deg in 13 days (8145 -> 8158): a flyby of Tellus at 0.03-0.05 AU. The orbits before (a=-5.7, e=1.14) and after (a=-2.5, e=1.30) the encounter both fit at 1.15 arcmin, same inclination 4.66 deg, so it is one object deflected by Tellus. Fitting the flyby with a point mass needs m_T = 1.29e-4 m_S and still only reaches 4.6 arcmin; C4 (and C7 over 3000 days) prefer m_T = 3-4e-6 (misfit rises 5x at 1e-5 and 20x at 3e-5). These two numbers disagree by a factor ~40. So the Tellus gravity law at short range (or the flyby model) is incomplete. This is an open failure of my laws.
4. *Why a model without Tellus's mass fails.* Near-Earth bodies (C4, C7) need m_T: fits of a 2500-day window of C4 give median 13.9 arcmin with m_T = 0, 4.1 with 3e-6; 23 with 1e-5.

## 2. Laws (equations and constants)

**Stage-1 laws (unchanged).** Newtonian inverse-square gravity between S, Tellus and P2-P5 in an inertial frame fixed against the stars; P1 and M as test particles (M around Tellus with S's tide); ecliptic frame tilted by eps = 21.702 deg; constants exactly as in ANSWER_1.md / `predict.py` (n_Tellus = 0.017202503 rad/day etc.). In the units used, G m_S = n_T^2 with Tellus's semi-major axis = 1 and time in days.

**New law for faint bodies.** Every faint body C_i is a massless test particle with heliocentric position r (ecliptic frame):

    d2r/dt2 = -(1 - beta_i) G m_S r/|r|^3  +  sum_{j in {T,P2,P3,P4,P5}} G m_j [ (r_j - r)/|r_j - r|^3  -  r_j/|r_j|^3 ]

with G m_S = n_T^2 = 2.959e-4 AU^3/day^2 (a_T = 1), m_T = 3.0e-6 m_S (from C4/C7; range 2.9-4.2e-6; not determined by S or P1-P5), m_P2 = 0 (|m_P2| < 4e-7 from stage 1), m_P3 = 1.5146e-4, m_P4 = 1.0003e-3, m_P5 = 2.491e-4 (stage 1, each +-5e-7 to 7e-7), r_j the stage-1 heliocentric positions. The second term includes the indirect term. P1 and M are not used as perturbers (their masses are not known; P1 is not needed to reach noise level in the fits above). Distances are relative (unit = Tellus's orbit radius); the direction seen from Tellus is that of r - r_T, rotated into (lon, lat) as in stage 1.

Constants of each body (osculating heliocentric state at the epoch; x,y,z in AU, vx,vy,vz in AU/day, ecliptic frame; 1-sigma from the least-squares covariance of the fitted window, which does not include the model errors described above; beta uncertainties are formal values <1e-4, the real spread between apparition pairs is 0.001-0.006):

| body | epoch (day) | x | y | z | vx | vy | vz | beta | a (AU), e, i (deg) |
|---|---|---|---|---|---|---|---|---|---|
| C1 | 8245 | -1.208297 | 1.023362 | 0.056640 | -0.02213722 | 6.2267e-05 | -4.8047e-04 | 0 (fixed) | -2.53, 1.298, 4.81 (hyperbolic; post-flyby) |
| C2 | 33137 | -2.595340 | 0.040428 | -0.401879 | 0.00139534 | -0.01170859 | 0.00451317 | 0.0007 | 4.496, 0.446, 21.64 |
| C3 | 5114 | -0.342603 | -1.709052 | 0.064346 | 0.01867796 | -0.00162304 | 0.00311560 | -0.254 +- 0.003 | 5.779, 0.702, 9.94 |
| C4 | 36071 | -1.716957 | 1.090242 | -0.016096 | -0.00604903 | -0.00707510 | 0.00013679 | 0.00001 | 1.448, 0.424, 1.02 |
| C5 | 30660 | 1.080160 | -2.451025 | -0.380833 | 0.01032638 | 0.00547821 | 0.00383790 | -0.0025 | 4.367, 0.394, 19.24 |
| C6 | 34685 | -0.660647 | -2.680939 | -0.751915 | 0.01126479 | -0.00246014 | 0.00119642 | 0 (fixed) | 4.085, 0.304, 16.12 |
| C7 | 35266 | -0.657729 | -0.706871 | 0.013743 | 0.01364401 | -0.01255448 | 0.00567111 | 0.0008 | 1.250, 0.228, 17.03 |
| C8 | 36071 | -0.066821 | 1.311636 | 0.283772 | -0.01427116 | -0.00842764 | 0.00161009 | +0.2971 | 6.431, 0.834, 17.24 |
| C9 | 20396 | 0.382853 | 2.109661 | -0.549265 | -0.01220029 | -0.00602091 | 0.00286280 | +0.0816 | 5.198, 0.764, 14.90 |
| C10 | 17169 | -1.411750 | 0.542150 | 0.281409 | 0.00695606 | -0.01768074 | -0.00021904 | -0.0256 | 9.043, 0.909, 13.61 |
| C11 | 35630 | -2.249923 | 0.065563 | 0.058244 | -0.00757115 | -0.01220356 | -0.00024848 | -0.150 (fixed from many apparitions) | 3.546, 0.593, 2.49 |
| C12 | 34661 | -1.844515 | -1.219802 | -0.092822 | 0.00979865 | -0.01370615 | 0.00472515 | -0.2962 | 9.491, 0.767, 15.76 |

(a, e, i are osculating at the epoch with mu = (1-beta) G m_S; 1-sigma on the state vector components are of order 1e-4 to 1e-3 AU in position and 1e-7 to 1e-5 AU/day in velocity for C2, C4, C5, C8, C9, C11, C12, and 1e-3 AU / 1e-5 AU/day for C3, C7; for C1 2e-3 AU.)

The orbits are integrated numerically (DOP853, rtol 1e-10) in `predict_2.py`, exactly as in the fits; stage-1 bodies by the stage-1 equations.

**What I could not determine.** (a) Whether beta is an inertia/gravity-coupling ratio or a radial non-gravitational force (the data cannot tell; I assumed only the term on S is scaled; planets' pull on faint bodies was tested only on C11, scale 1.005); (b) the strength of Tellus's gravity at short range (C1 versus C4/C7); (c) masses of P1, P2 and M (they affect only encounters; stage 1 cannot see them); (d) the true orbits of C6 (only one apparition fits), C3 (last data day 14297), C10, C1.

## 3. Predictions

`predictions_2.csv` was produced by `python predict_2.py` using exactly the constants above and the stage-1 constants in `predict.py`.
Expected quality: stage-1 bodies as in stage 1 (a few arcmin). C2, C4, C5, C8, C9, C12: I expect errors of order arcmin to ~1 deg unless an unmodelled encounter happens. C11: last orbit uses beta=-0.150, expect up to a few degrees. C3 (extrapolated 22,000-29,000 days beyond its last data), C6, C10, C7: errors of many degrees possible. C1: after a hyperbolic flyby it is 320-400 AU away; its direction should be good to ~0.5 deg.

## 4. Self-assessment

| statement | confidence |
|---|---|
| Stage-1 laws (Newtonian N-body, eps = 21.702 deg) hold for S, P1-P5, M and are unchanged | 95 % |
| Faint bodies are massless test particles obeying inverse-square attraction by S plus the planets (perturbations included) | 80 % |
| Each faint body has a constant factor (1-beta_i) on S's attraction, with beta_i as tabulated (C11, C12, C8, C9, C3 clearly non-zero; C2, C4, C5, C7 zero) | 70 % for existence of non-zero beta; 50 % that "constant per body, only on S's term" is the exact form |
| m_T = 3e-6 m_S | 40 % (consistent between C4 and C7, in conflict with C1) |
| Orbits of C2, C4, C5, C8, C9, C12 will be right to < 1 deg in years 101-120 | 60 % each |
| C1, C3, C6, C7, C10, C11 predictions right to < 5 deg | 35 % |

What would prove me wrong: (1) A faint body whose recent orbit shows an apparition in years 101-120 at a time or place that differs by more than a few days / a few degrees from my prediction (periods: C2 about 9.5 yr, C5 9.1 yr, C12 25.6 yr, C8 ~19 yr, C9 ~12 yr). (2) Another close Tellus flyby whose deflection is inconsistent with m_T = 3e-6 (as C1). (3) beta values that change from apparition to apparition beyond 0.005 (then beta is not a constant of the body). (4) A planetary encounter (P3-P5) whose observed deflection disagrees with the stage-1 masses.
