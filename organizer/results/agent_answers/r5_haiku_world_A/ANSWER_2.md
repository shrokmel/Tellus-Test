# ANSWER_2

## 1. Stage 1 laws tested against all bodies

**Holds (within noise):** S, P1, P4, P5 (and P2 apart from a small drift). Kepler orbits about the Sun with the common GM hold to about 1 arcmin per component; per-window residuals 0.8-1.2 arcmin. P3 is held to 2-4 arcmin with an empirical harmonic correction (stage 1 residual 8 arcmin). Moon: stage 1 harmonic model, 18 arcmin in-sample.

**Holds for individual arcs of faint bodies, but not for whole records:**
- C2 arcs 1-3 (6456-6673, 6783-6970, 9507-10049) each fit a heliocentric Kepler orbit at noise level (0.79-0.82 arcmin rms), with a^3 n^2 / n_E^2 = 0.99-1.00 (Law 1 holds per arc). Their inclinations (0.3786 rad) and nodes (about 3.136 rad) match Earth's plane (0.3788 rad, 3.142 rad).
- But a and e change between arcs: arc 1 a=7.24 AU, e=0.78; arc 2 a=4.32, e=0.62; arc 3 a=4.14, e=0.61. Elements that jump between arcs are not consistent with a single fixed Kepler orbit.

**Fails for whole records (rms residual in arcmin; noise is about 0.8 arcmin per component):**

| body | sightings | heliocentric Kepler, common GM | heliocentric, free a | geocentric Kepler |
|---|---|---|---|---|
| C1 | 205 | 97 | 76 | 963 |
| C2 | 2873 | 1478 | 1602 | 1562 |
| C3 | 989 | 242 | - | - |
| C4 | 14425 | 1938 | - | - |
| C5 | 3271 | 1341 | - | - |
| C6 | 434 | 368 | - | - |
| C7 | 21900 | 2141 | 2124 | 1884 |
| C8 | 2009 | 323 | - | - |
| C9 | 2502 | 853 | - | - |
| C10 | 762 | 310 | - | - |
| C11 | 6289 | 1254 | - | - |
| C12 | 885 | 634 | - | - |

Search check: the same random-start search recovers P2 (1.4 arcmin) from scratch, so the failures are not a search failure. Other stage 1 fitting failures (for example P1 with the wrong period grid) were solved only when I searched a wider period range.

**Conclusion of the test:** the stage 1 laws (single Kepler orbit, common GM) do not describe the faint bodies over the whole record. The faint bodies' orbits change between arcs, or they are perturbed (close approaches to planets, or a force not in stage 1). I could not finish a per-arc analysis for all arcs; the per-arc sweep was too slow to complete.

## 2. Laws (equations and constants)

Laws 1-5 of stage 1 are retained for S, P1-P5 and M (see ANSWER_1.md; constants unchanged).

**Law F (faint bodies, provisional).** Each faint body C_k moves on a heliocentric Kepler orbit with the common GM of Law 1:
  r_C(t) = Kepler(a_k, e_k, i_k, Omega_k, omega_k, M0_k, n_k; t), with a_k = (GM_S/n_k^2)^(1/3), GM_S = n_E^2 = 2.9598e-4 AU^3/d^2.
  Geocentric direction: (r_C - r_E)/|r_C - r_E|, as in Law 3 of stage 1.
  Whole-record fitted constants (in predict_2.py):

| body | e | i [rad] | Omega [rad] | omega [rad] | M0 [rad] | n [rad/d] | a [AU] |
|---|---|---|---|---|---|---|---|
| C1 | 0.6968 | 0.3310 | -3.1020 | 3.9773 | 19.8203 | 0.0037193 | 2.776 |
| C2 | 0.4388 | 0.5339 | 3.2052 | 3.4620 | 5.5963 | 0.0090014 | 1.540 |
| C3 | 0.3832 | 0.4634 | 3.4768 | 1.5172 | 3.5270 | 0.0041932 | 2.563 |
| C4 | 0.9500 | 0.3953 | 3.1386 | 2.8131 | 0.6266 | 0.0223294 | 0.840 |
| C5 | 0.7258 | 0.3401 | 3.2547 | 6.1019 | -3.9772 | 0.0020501 | 4.129 |
| C6 | 0.8441 | 0.1424 | 3.6044 | 4.5855 | 10.9215 | 0.0018915 | 4.357 |
| C7 | 0.9500 | 2.0798 | 6.8990 | 2.8415 | 5.0220 | 0.0771525 | 0.368 |
| C8 | 0.5788 | 0.3259 | 2.0965 | 6.7347 | 3.1794 | 0.0035630 | 2.857 |
| C9 | 0.2465 | -0.6175 | 6.0793 | 3.0701 | 8.0006 | 0.0069763 | 1.825 |
| C10 | 0.5905 | 0.5114 | 2.6889 | 1.6883 | 5.9475 | 0.0051146 | 2.245 |
| C11 | 0.5158 | -0.3439 | 6.1878 | 1.6733 | 7.2004 | 0.0032083 | 3.063 |
| C12 | 0.7385 | 0.6075 | 3.4692 | 0.7355 | 7.7285 | 0.0006825 | 8.597 |

Uncertainties: the rms residuals in section 1 are the honest error scale. These constants describe the whole-record average and do not describe any arc of any faint body well (for C2, per-arc a and e differ by up to 0.16 and 0.2 from the whole-record values).

**Note on the stage 1 laws for the rest:** S, P1-P5 and M use the stage 1 equations and constants exactly as in `predict.py`, which `predict_2.py` imports.

## 3. Predictions

`predictions_2.csv` (all 4617 rows, all 19 bodies, years 101-120) is filled by running `predict_2.py`. S, P1-P5 and M use the stage 1 predictor; C1-C12 use the provisional Law F. The predictions for C1-C12 are extrapolations from whole-record averages. They are not expected to match the faint bodies, and I would not rely on them.

## 4. Confidence and what would prove me wrong

- **Stage 1 laws for S, P1, P2, P4, P5 (Kepler, common GM), about 90%.** They hold to noise over 100 years. Proof of error: a residual above 3 arcmin in a new window for any of these bodies.
- **Stage 1 P3 and Moon harmonic models, about 50-60%.** Both are empirical stand-ins for perturbations. Proof of error: out-of-sample residuals above 30 arcmin.
- **Law F (faint bodies, whole-record Kepler), about 5%.** Whole-record residuals are 97 to 2141 arcmin, and the per-arc elements differ between arcs. It is not a law. Proof that would make it right: a single heliocentric orbit fitting all sightings of a body to within about 1-2 arcmin rms.
- **The thing that would change the picture:** a physical law for the faint bodies that stage 1 does not have, for example a non-gravitational force or a much larger GM for light bodies, or perturbations by planets at close approaches. A test would be to compute each body's close approaches to the stage 1 planet orbits and see whether the element jumps in C2 (and other faint bodies) occur at those times.

Status: stage 2 is not complete. The per-arc analysis for all bodies is unfinished, and I did not find the law that explains the faint bodies.
