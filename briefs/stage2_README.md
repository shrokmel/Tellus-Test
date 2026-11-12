# Stage 2: the faint bodies

Your answer to stage 1 has been saved and will not be changed. The observers on Tellus have now released their
records of the faint bodies in the sky, over the same 100 years and with the same precision, gaps and format.

## New data (`observations_2.csv`)
- `C1`–`C12`: twelve faint bodies. A faint body is recorded only on nights when it is bright enough to be seen, so
  most are seen only for some weeks or months at a time.
- `observations.csv` from stage 1 is still here.

## Your tasks
1. **Test your stage 1 laws** against all the bodies. Say plainly where they hold and where they fail, with
   numbers.
2. **Give the laws as equations.** State the simplest set of laws that explains every body, old and new, at
   once, as explicit equations with a numerical value (and uncertainty) for every constant. If they are the same
   as in stage 1, say so; if they differ, give the new equations and say what forced the change. Words alone do
   not count.
3. **Predict** the direction of every body on the days listed in `predictions_2.csv` (years 101–120), using a
   script `predict_2.py` that uses only the equations and constants you stated in task 2 and the observations in
   this folder. We will run it. Predictions that do not come from your stated equations do not count.
4. **Judge yourself**: how confident are you in each law (0–100%)? What observation would prove you wrong?

## Rules
Same as stage 1. Write your answer to `ANSWER_2.md` and complete `predictions_2.csv` by running `predict_2.py`.
