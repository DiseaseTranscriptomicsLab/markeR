# Compute Statistical Tests for Variable Associations with a Target Variable

Performs statistical tests to assess the relationship between predictor
variables and a target variable, selecting appropriate methods based on
variable types. Returns a list of data frames containing metric values
and p-values.

## Usage

``` r
compute_stat_tests(
  df,
  target_var,
  cols = NULL,
  numeric = "pearson",
  categorical_bin = "t.test",
  categorical_multi = "anova",
  p.adjust.method = "BH"
)
```

## Arguments

- df:

  A data frame containing the target variable and predictors.

- target_var:

  A string specifying the dependent variable.

- cols:

  Optional. A character vector of predictor variables. If `NULL`, all
  variables except `target_var` are used.

- numeric:

  The correlation method for numeric predictors. Options: `"pearson"`
  (default), `"spearman"`, `"kendall"`.

- categorical_bin:

  The statistical test for binary categorical variables. Options:
  `"t.test"` (default) or `"wilcoxon"`.

- categorical_multi:

  The statistical test for multi-level categorical variables. Options:
  `"anova"` (default) or `"kruskal-wallis"`.

- p.adjust.method:

  Character string specifying the method to use for multiple testing
  correction. Must be one of `"BH"` (Benjamini-Hochberg, default),
  `"holm"`, `"hommel"`, `"bonferroni"`, `"BY"` (Benjamini-Yekutieli),
  `"fdr"`, or `"none"`. Passed to
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html).

## Value

A named list (one entry per variable being analysed) where each element
is a data frame with:

- **Metric**: The test statistic (correlation coefficient, t-statistic,
  ANOVA F-value, etc.).

- **p-value**: The significance value of the test.

- For **Categorical Multi**, multiple rows are included for pairwise
  comparisons (Tukey HSD results).

## Details

### **Variable Classification**

- **Numeric**: Continuous numeric or integer variables with more than 10
  unique values.

- **Categorical Bin**: Binary categorical variables (factors,
  characters, or integers with exactly 2 unique values).

- **Categorical Multi**: Categorical variables with more than 2 unique
  values (up to 10 levels recommended). A warning is issued for
  categorical variables with more than 10 unique values.

### **Statistical Tests Applied**

- **Numeric Predictors**: Pearson, Spearman, or Kendall correlation.

- **Categorical Bin Predictors**: T-test or Wilcoxon rank-sum test.

- **Categorical Multi Predictors**: ANOVA (default) or Kruskal-Wallis
  test. If ANOVA is used, Tukey's HSD post-hoc test is performed for
  multiple comparisons.

The function automatically detects variable types and applies the
appropriate test. If a categorical variable has more than 10 unique
levels, a warning is issued. If an invalid statistical test is
requested, the function stops with an error message.

## Examples

``` r
if (FALSE) { # \dontrun{
df <- data.frame(
  score = c(80, 85, 90, 95, 100),
  age = c(25, 30, 35, 40, 45),
  gender = c("Male", "Female", "Male", "Female", "Male"),
  group = factor(c("A", "B", "A", "B", "C"))
)

results <- compute_stat_tests(df, target_var = "score")
print(results)
} # }
```
