# Generate All Possible Unique Contrasts Between Groups

This function creates statistical contrasts between levels of a
categorical variable. Users can choose the level of complexity:

- `"simple"`: Pairwise comparisons (e.g., A - B).

- `"medium"`: Pairwise comparisons plus comparisons against the mean of
  other groups.

- `"extensive"`: All possible groupwise contrasts, ensuring balance in
  the number of terms on each side.

## Usage

``` r
generate_all_contrasts(levels, mode = "simple")
```

## Arguments

- levels:

  A character vector of unique group levels.

- mode:

  A string specifying the level of detail for contrasts. Options are
  `"simple"` (pairwise only), `"medium"` (pairwise + vs. mean of
  others), or `"extensive"` (all possible balanced groupwise contrasts).
  Default is `"extensive"`.

## Value

A character vector of unique contrast expressions.

## Examples

``` r
if (FALSE) { # \dontrun{
levels <- c("A", "B", "C", "D")
generate_all_contrasts(levels, mode = "simple")    # Pairwise only
generate_all_contrasts(levels, mode = "medium")    # Pairwise + mean comparisons
generate_all_contrasts(levels, mode = "extensive") # All balanced contrasts
} # }
```
