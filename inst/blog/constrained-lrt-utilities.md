# A Small Utility Package for Constrained LRT Reference Laws

Sometimes the hard part of a constrained likelihood ratio test is not fitting the model. It is what comes after: figuring out how to work with the reference distribution once it is no longer a single chi-square law.

That is where `mixchisq` comes in.

The package is intentionally small. It gives you four basic tools for finite mixtures of chi-square distributions:

- `dmixchisq()`
- `pmixchisq()`
- `qmixchisq()`
- `rmixchisq()`

and a quick visual helper, `plot_mixchisq()`.

## Why this is useful in practice

In ordinary likelihood ratio testing, you often compare your test statistic to a single chi-square distribution. But in constrained settings, that can change. The limiting law may be a mixture instead. In applied work, that usually means one very practical problem: you need to evaluate probabilities and critical values for a distribution that base R does not provide directly.

A common example is when a constrained LRT produces a Self-and-Liang-style limiting law. In that case, you may end up with something like:

- `0.5 * chi^2_0 + 0.5 * chi^2_1`, or
- `0.5 * chi^2_2 + 0.5 * chi^2_3`.

At that point, you do not necessarily need a long theoretical derivation in your codebase. You just need a dependable utility that can answer questions like:

- What is the CDF at this observed statistic?
- What is the 95th percentile?
- How do I simulate from this law for a calibration plot?
- Can I replace my hand-written mixture simulation code with one clean function call?

That is the role of `mixchisq`.

## A few places it helps

### 1. Constrained LRT workflows

If your asymptotic reference law is already known, you can use the package immediately for p-values and critical values.

### 2. Simulation studies

If you are checking empirical LRT behavior against a mixture reference law, `rmixchisq()` and `pmixchisq()` simplify the code considerably.

### 3. Reproducible applied analysis

A named function is easier to audit than an ad hoc block of simulation code embedded inside an analysis script.

### 4. Teaching and communication

Plots and direct function calls make it easier to explain what reference law you are using without dragging readers through implementation details.

## Example

Instead of writing a manual latent-mixture simulation every time, you can do:

```r
rmixchisq(10000, df = c(2, 3), w = c(0.5, 0.5))
```

And if you need a critical value:

```r
qmixchisq(0.95, df = c(2, 3), w = c(0.5, 0.5))
```

That is the whole idea. Keep the package narrow, useful, and easy to drop into applied work.

## Final thought

`mixchisq` is not trying to replace asymptotic theory. It is a convenience layer for the moment after the theory has already told you what limiting law to use. Once you know the mixture, the package helps you work with it cleanly.
