# ANSWER 2 — The faint bodies

Short version:
- **Gravity is still Newton's inverse-square law, and the bright bodies still follow it exactly.** Two things change.
- **Tellus has a measurable mass**, GM_T/GM_S = (3.0 ± 0.5) × 10⁻⁶. In stage 1 I had set it to zero.
- **Each faint body responds to gravity with its own strength:** its acceleration is the Newtonian one multiplied by a constant κᵢ = 1 − βᵢ, with βᵢ between −0.30 and +0.30. This factor applies to the pull of S and of every other mass alike, and it is the only new constant per body.
- **βᵢ is fixed for decades but can jump at a close encounter** with Tellus or a P (seven such events found). The predictions use each body's value after its last encounter.

Scripts: `comet.py`, `ng.py`, `nglink2.py`, `backlinkK.py`, `c1fit.py`, `c7fit.py`, `refit_recent.py` (fits); `segtab.py`, `perapp.py`, `linktest.py`, `tseg.py`, `check_2.py`, `unc_2.py` (tests); `make_predict_2.py` (writes the constants into `predict_2.py`); `predict_2.py` (predictions).

---

## 1. Testing the stage-1 laws on all bodies

**Bright bodies (S, P1–P5, M).** The stage-1 laws hold. With the stage-2 constants (Tellus's mass included, section 2) the residuals are S 1.00′, P1 1.01′, P2 1.00′, P3 1.00′, P4 1.01′, P5 1.00′ and M 1.00′ rms per angle (`check_2.py`).

**Faint bodies: stage-1 laws as they stand (point masses S, P1–P5, Tellus mass 0) fail.** Each C was treated as a test particle in the stage-1 field. The table gives the fit to one apparition, using only a few months of data (`cseg.py`, `segtab.py`), and then to the whole record.

| body | gravity-only rms, best single apparition | gravity-only, whole record |
|---|---|---|
| C1 | 12.2′ | 12.2′ (single passage) |
| C2 | 1.08′ | fails (tens of arcmin to degrees) |
| C3 | 22.1′ | fails |
| C4 | 1.02′ | fails |
| C5 | 18.4′ | fails |
| C6 | 1.08′ | fails |
| C7 | 20′ (89′ for days 616–1612) | fails |
| C8 | > 100′ | fails |
| C9 | 19.9′ (36.5′ for days 2207–2705) | fails |
| C10 | 4.6′ | fails |
| C11 | 29.3′ | fails |
| C12 | 26.5′ | fails |

So the laws fail badly for 9 of the 12 bodies even over a single apparition of a few months. For those three that fit at about 1′ (C2, C4, C6), they still fail when apparitions years apart are joined.

**Diagnosis.**

1. **One extra radial 1/r² term fixes every apparition.** Within one apparition, adding a radial acceleration −β GM_S r̂/r² brings every apparition of every body (except those containing a close encounter) to 0.94–1.05′. Other radial dependences were tested: 1/r, 1/r³, a constant, and the comet outgassing law g(r). All are clearly worse; for example C8 gives 16′ for 1/r, 1.04′ for 1/r² and 13′ for 1/r³. Transverse and normal terms are not needed.
2. **β is a constant of each body.** It is the same in every apparition. C9 gives 0.0816–0.0820 in all 8 of its apparitions over 100 years; C3 gives −0.2303 to −0.2307; C7 gives 0.2173–0.2176 in three separate windows.
3. **The factor must also scale the planets' pull.** With β acting on S's pull only, two apparitions of C9 joined across one 12.4-year orbit give 2.2′. Scaling the direct pull of *every* mass by the same κ gives 0.99′, and it predicts the next, unfitted apparition to 1.9′. A model with separate factors for S and for the planets does not settle down (the planet factor wanders).
4. **Tellus's mass is needed.** Over many orbits the joined fits still drifted (C9: 1.0′ → 17′ over 8 orbits) until Tellus's mass was added. C7 orbits near Tellus. In its windows with β = 0 the residuals were ~1.0′ except in the windows with a close pass by Tellus (0.04–0.26 L), which gave 1.3′–6.5′. With GM_T free those windows fit to 0.99–1.02′ with GM_T/GM_S = 3.33 (±0.04), 2.93 (±0.12) and 2.85 (±0.30) × 10⁻⁶. One window gives 7.0 × 10⁻⁶; β probably also changed during that pass.

   The stage-1 planets accept, and slightly prefer, a Tellus mass: refitting stage 1 with GM_T/GM_S = 3.0 × 10⁻⁶ fixed gives 1.00227′ (zero mass: 1.00233′). With the mass free the planets give (3.5 ± 0.3) × 10⁻⁶, but that fit only converges by giving P1 a negative mass, so I did not adopt it. My stage-1 value (zero) was a fit that had not converged.

   With GM_T/GM_S = 3.0 × 10⁻⁶ included, C9's entire 100-year record (2,502 sightings, 8 orbits) fits with one constant β to 1.015′. C12's (885 sightings) fits to 0.996′ and C10's (762) to 0.985′.
5. **β changes at close encounters.** Some bodies' β changes abruptly, always at a close encounter. Within each stretch between encounters β is constant to ±10⁻⁴ or better:

| body | event | β before | β after |
|---|---|---|---|
| C1 | passes 0.03 L from Tellus, day ≈ 8151 | −0.263 ± 0.005 (days 8012–8147, 1.06′) | −0.0139 ± 0.0002 (days 8155–8331, 1.02′) |
| C5 | passes 0.054 L from P4, day ≈ 19270 | −0.1718 (6 apparitions, days 1930–17700, 1.08′) | −0.0010 ± 0.0001 (days 20606–34296, 1.05′) |
| C7 | 0.16 L from Tellus, day ≈ 4600 | 0.2176 | 0.105 |
| C7 | 0.07 L from Tellus, day ≈ 11000 | 0.105 | 0.0002 |
| C7 | (around day 32000) | 0.0002 | −0.021 |
| C7 | 0.07 L from Tellus, day 35777 | −0.021 | −0.0001 ± 0.0000 |
| C3 | day ≈ 14100 | −0.2305 (days 380–9834) | −0.063; not seen after day 14297 |
| C2 | orbit change between days 16020 and 19035 (perihelion 1.45 → 2.5) | ≈ +0.003 | 0.0001 ± 0.0001 |
| C11 | days 33330–33390 | ≈ −0.15 | −0.1155 |

Other tests that failed:
- light-time with a free speed of light (1/c came out negative);
- a hidden body on Tellus's orbit, opposite Tellus;
- a cut-off of the extra force beyond r = 3 (it contradicts C9's observed 4,500-day period);
- a static distant tidal field, which fits but needs a strength the planets clearly rule out.

**Visibility.** A faint body is recorded only when it is less than 3.0 L from S. This is a brightness threshold, not a law of motion. For example C2's r always lies within 1.45–3.00 and C7's within 0.96–3.0.

---

## 2. The laws as equations

Notation and units are as in stage 1. L is the length unit (Tellus's mean distance ≈ 0.9999 L). Time is in days from day 0. **r** is a heliocentric ecliptic position.

**Law 1 — frame (unchanged).** **e** = R(ε) **u**, with ε = 21.7015° ± 0.0010°. Directions are geometric, from Tellus.

**Law 2 — Newtonian gravity for S, Tellus and P1–P5 (same equation; the constants are refitted because m_T ≠ 0).**

d²**r**ᵢ/dt² = −G(M_S + mᵢ) **r**ᵢ/|**r**ᵢ|³ + Σ_{j≠i} G mⱼ [ (**r**ⱼ − **r**ᵢ)/|**r**ⱼ − **r**ᵢ|³ − **r**ⱼ/|**r**ⱼ|³ ],  i, j ∈ {T, P1, …, P5}

| constant | value (≈ 2σ formal) |
|---|---|
| G M_S | 2.958466 × 10⁻⁴ L³ d⁻² (defines L) |
| m_T/M_S (Tellus + M) | **3.0 × 10⁻⁶ ± 0.5 × 10⁻⁶** (new; from C7's passes by Tellus and C9; fixed at 3.0 × 10⁻⁶ in the fit) |
| m_P1/M_S | (0.57 ± 0.17) × 10⁻⁶ |
| m_P2/M_S | (0.94 ± 0.24) × 10⁻⁶ |
| m_P3/M_S | (1.4974 ± 0.0010) × 10⁻⁴ |
| m_P4/M_S | (1.00077 ± 0.00011) × 10⁻³ |
| m_P5/M_S | (2.4884 ± 0.0014) × 10⁻⁴ |

Starting state at t = 0 (L and L/day):

| body | x | y | z | vx | vy | vz |
|---|---|---|---|---|---|---|
| Tellus | -0.394274928 ± 7.5e-06 | -0.938258670 ± 6.4e-06 | -0.000024789 ± 3.7e-06 | 0.01551748490 ± 1.3e-07 | -0.00668736080 ± 1.0e-07 | -0.00000017848 ± 6.1e-08 |
| P1 | -0.312554735 ± 1.1e-05 | -0.301676255 ± 1.6e-05 | -0.000003525 ± 5.9e-06 | 0.01747714234 ± 6.8e-07 | -0.01997912525 ± 8.0e-07 | -0.00000019404 ± 3.6e-07 |
| P2 | -0.714792200 ± 1.4e-05 | 1.633058239 ± 1.4e-05 | 0.000042255 ± 8.0e-06 | -0.01145284526 ± 8.8e-08 | -0.00516815528 ± 8.5e-08 | -0.00000009348 ± 5.6e-08 |
| P3 | -3.056442389 ± 2.6e-05 | 1.364970172 ± 3.0e-05 | 0.053062029 ± 1.7e-05 | -0.00411729757 ± 7.5e-08 | -0.00837383641 ± 6.5e-08 | 0.00019967811 ± 4.8e-08 |
| P4 | 6.460515141 ± 5.0e-05 | -0.028485271 ± 6.5e-05 | -0.032006542 ± 3.4e-05 | 0.00021088049 ± 6.0e-08 | 0.00662642870 ± 5.3e-08 | -0.00011165154 ± 3.5e-08 |
| P5 | 3.323331349 ± 1.0e-04 | 11.027411654 ± 8.7e-05 | 0.000175215 ± 6.1e-05 | -0.00490180861 ± 3.7e-08 | 0.00126579410 ± 4.4e-08 | 0.00017668013 ± 2.7e-08 |

Full-precision values are in `predict_2.py`. Fit: 1.0023′ rms over all S and P sightings. Integration: RK4, step 0.25 d.

**Law 3 — M about Tellus (unchanged).** The equation and constants are exactly as in stage 1. Directions of M fix only the month and the shape of the orbit, so Tellus's new mass does not alter them.

**Law 4 — faint bodies (new).** Each faint body Cᵢ is a massless body. Its acceleration is the Newtonian field of all the masses, multiplied by its own constant κᵢ:

d²**r**/dt² = κᵢ [ −G M_S **r**/|**r**|³ + Σ_{j∈{T,P1..P5}} G mⱼ (**r**ⱼ − **r**)/|**r**ⱼ − **r**|³ ] − Σ_{j∈{T,P1..P5}} G mⱼ **r**ⱼ/|**r**ⱼ|³,  κᵢ = 1 − βᵢ

The last sum is the usual heliocentric frame term (S's own acceleration) and is not scaled. For the bright bodies κ = 1: their Kepler's-third-law ratios agree to 10⁻⁵ and the N-body fit is at the noise level.

βᵢ is constant between close encounters and changes abruptly at a close encounter (distance ≲ 0.2 L from Tellus or a P). Section 1 lists the observed changes. I have no law for the size of the change; that is the open part of this law.

Constants for the latest stretch of each body (used for prediction). State at epoch t₀ is heliocentric ecliptic, L and L/day. β is ± 2σ formal.

| body | t₀ (day) | x | y | z | vx | vy | vz | β | arc fitted (days) | rms |
|---|---|---|---|---|---|---|---|---|---|---|
| C1 | 8245 | -1.18865939 | 1.02254762 | 0.05556596 | -0.0222390629 | 0.0000377491 | -0.0004633069 | -0.013892 ± 0.000446 | 8155–8331 | 1.016′ |
| C2 | 33237 | -2.23741746 | -1.09902997 | 0.07082077 | 0.0057970908 | -0.0107159842 | 0.0048029618 | +0.000063 ± 0.000152 | 19035–33490 | 1.018′ |
| C3 | 14212 | 1.91432326 | -1.08102914 | 0.37188786 | 0.0142140547 | 0.0084591245 | 0.0014507675 | -0.063463 ± 0.010605 | 14126–14297 | 0.980′ |
| C4 | 35935 | -0.38758609 | 1.52539240 | -0.02628202 | -0.0128842902 | 0.0023612778 | -0.0000175975 | +0.000082 ± 0.000033 | 30454–36282 | 1.025′ |
| C5 | 34065 | 1.80356970 | -1.91982573 | -0.07283844 | 0.0085848464 | 0.0081397814 | 0.0041165293 | -0.001013 ± 0.000194 | 20606–34296 | 1.048′ |
| C6 | 34740 | -0.03204879 | -2.75977318 | -0.67080632 | 0.0115063444 | -0.0005364930 | 0.0017001735 | +0.035152 ± 0.006448 | 34575–34907 | 0.962′ |
| C7 | 36112 | 0.15288570 | 1.41215277 | -0.27735645 | -0.0128282457 | -0.0013435499 | -0.0024511990 | -0.000127 ± 0.000060 | 35781–36525 | 0.976′ |
| C8 | 36123 | -0.75647602 | 0.70298229 | 0.32045195 | -0.0112262666 | -0.0150297646 | -0.0004433230 | +0.297354 ± 0.000072 | 21803–36402 | 1.015′ |
| C9 | 33943 | -1.04275699 | 0.73891812 | -0.03053544 | -0.0077737852 | -0.0170569725 | 0.0046231408 | +0.081814 ± 0.000027 | 2207–34221 (all) | 1.015′ |
| C10 | 27114 | 0.03638771 | -0.85675561 | 0.06173673 | 0.0250384116 | -0.0043198738 | -0.0056078414 | -0.032576 ± 0.000236 | 7265–27301 (all) | 0.985′ |
| C11 | 35629 | -2.22082821 | 0.08372913 | 0.05806443 | -0.0075991051 | -0.0120619063 | -0.0002459922 | -0.115529 ± 0.001258 | 33390–35741 | 0.939′ |
| C12 | 34671 | -1.74383131 | -1.35504762 | -0.04548749 | 0.0104368107 | -0.0132485127 | 0.0047505871 | -0.297015 ± 0.000273 | 6909–34857 (all) | 0.996′ |

Full precision is in `predict_2.py`. Law 4 is integrated with DOP853 (relative tolerance 10⁻¹¹), with the Law-2 bodies tabulated every 0.25 d and Hermite-interpolated.

**What forced the changes.**
- **Tellus's mass:** bodies passing close to Tellus (C7 and C1) and the many-orbit drift of C9. The planets are consistent with it.
- **κᵢ:** single apparitions of 9 of the 12 faint bodies are off by 4′–170′ under pure gravity. No mass distribution, frame effect or light-time effect can do that while leaving the planets untouched. A per-body factor on the response to gravity fixes all apparitions with one number each.
- **β changes at encounters:** arcs on either side of a close encounter need different β, each fitting to the noise.

---

## 3. Predictions

`predict_2.py` fills `predictions_2.csv` (4,617 rows) in about 30 s. It uses:
- S, P1–P5: Law 2 with the stage-2 constants;
- M: Law 3;
- C1–C12: Law 4 with the constants above, integrated from each body's epoch t₀.

Expected accuracy in years 101–120. The formal 1σ values are from 12 random draws from each body's covariance (`unc_2.py`); the median is first, the worst case in brackets.

| body | expected accuracy |
|---|---|
| S, P1–P5, M | as in stage 1 (< 0.1′; P1 ≲ 1′) |
| C2, C4, C5, C8, C9, C10, C12 | 0.04–0.2′ (≤ 1.4′) |
| C1 | 1.3′ |
| C11 | 0.5′ (4.6′) |
| C7 | 0.8′ (53′ around close passes by Tellus) |
| C3 | ~95′, based on one 171-day window 60+ years ago |
| C6 | ~200′, based on one 332-day window near perihelion |

These are formal errors only. The real risk is a β change at a future close encounter, which Law 4 cannot foresee. C7 meets Tellus every few years and C4 passes inside Tellus's orbit, so those two are most exposed.

---

## 4. Self-assessment

| law | confidence | what would prove it wrong |
|---|---|---|
| Law 1 (frame) | 98% | as stage 1 |
| Law 2 (Newtonian N-body, now with m_T = 3.0 × 10⁻⁶ M_S) | 93% for the form; 75% for m_T to ±0.5 × 10⁻⁶ | any P or S departing from `predictions_2.csv` by ≫ 1′; a near-Tellus faint body deflected by a different amount than this mass predicts |
| Law 3 (M) | 90% | as stage 1 |
| Law 4a: faint body acceleration = κᵢ × Newtonian field, κᵢ = 1 − βᵢ, same factor for all masses | 75% | a faint body whose single apparition cannot be fitted with a constant κ; a close planetary encounter whose deflection needs a different factor for the planet than for S |
| Law 4b: βᵢ constant between close encounters | 70% | C9, C10, C12 (no recent encounter) drifting from prediction by more than a few arcmin; a β change with no encounter |
| Law 4c: β changes only at close encounters (≲ 0.2 L) | 50% | a change seen far from any encounter, or a close encounter with no change. The size of the jump is not explained, so this is the least certain part |

Weak points:
- The mechanism behind κ ≠ 1 is unknown. Physically it looks like a body-specific ratio of gravitational to inertial mass. Radiation pressure would give only β > 0, but β < 0 is observed (C3, C5, C10, C11, C12, C1).
- The rule for β changes is empirical and incomplete.
- C3 and C6 have poorly constrained orbits.
- The Tellus mass is known only to about 15%.
