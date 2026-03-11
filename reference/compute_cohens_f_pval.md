# Compute Cohen's f and p-value for a given model and predictor

This function calculates the Cohen's f effect size and the corresponding
p-value for a given linear model or ANOVA model based on the predictor
variable type (numeric or categorical).

## Usage

``` r
compute_cohens_f_pval(model, type)
```

## Arguments

- model:

  A linear model (`lm`) or ANOVA model (`aov`) fitted to the data.

- type:

  A string indicating whether the predictor is numeric or categorical.
  Options are "Numeric" or "Categorical".

## Value

A named vector with two elements:

- `Cohen_f`: The Cohen's f effect size value.

- `P_Value`: The p-value from the statistical test.

## Details

Cohen's f effect size is computed from the eta-squared (\\\eta^2\\)
value. For numeric predictors (continuous variables), the p-value is
obtained from the t-test in `summary(lm(...))`. For categorical
predictors (binary or multi-level), the p-value is obtained from the
F-test in `anova(lm(...))`.
