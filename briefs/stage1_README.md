# Observations of the sky from Tellus

You live on Tellus. This folder holds 100 years of nightly observations of the sky made from Tellus. Your task is
to find the laws that govern the motions you see, write them as equations, and use them to predict.

## The data (`observations.csv`)
- One row per sighting: `day`, `body`, `lon`, `lat`.
- `day` counts days since observations began (one day = one rotation of Tellus).
- `lon`, `lat` give the body's direction on the sky in degrees, in a fixed frame against the distant stars. The
  effect of Tellus's daily rotation has already been removed.
- Measurement error is about 1 arcminute in each angle.
- A body cannot be seen when it is within 15° of `S` on the sky, or when the sky is cloudy (about 30% of nights).
  `S` itself is measured by day, and about 30% of days are cloudy. So there are gaps.
- Bodies: `S` (the bright body that makes daylight), `P1`–`P5` (wandering bodies), `M` (a large bright body).
- Only directions are measured. No distances, sizes or brightnesses are given.

## Your tasks
1. **Describe** the motions you find, and how they change over the 100 years.
2. **Give the laws as equations.** State the simplest set of laws that explains the motions of all the bodies at
   once, as explicit equations with a numerical value (and uncertainty) for every constant. Words alone do not
   count.
3. **Predict** the direction of every body on the days listed in `predictions.csv` (years 101–120). Your
   predictions must be produced by a script, `predict.py`, that uses only the equations and constants you stated
   in task 2 and the observations in this folder. We will run it. Predictions that do not come from your stated
   equations do not count.
4. **Judge yourself**: how confident are you in each law (0–100%)? What observation would prove you wrong?

## Rules
- Work only inside this folder. Python with numpy and scipy is installed. Write whatever code you need.
- No internet access, and no files outside this folder.
- Write your answer to `ANSWER_1.md` (sections 1–4), and complete `predictions.csv` by running `predict.py`.
