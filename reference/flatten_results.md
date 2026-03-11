# Flatten a Nested List of Results into a Data Frame

Converts a nested list (where the first level is a method, the second
level is a gene signature, and the third level is a data frame) into a
single data frame. Additional columns for method and signature are added
to the data frame.

## Usage

``` r
flatten_results(nested_list)
```

## Arguments

- nested_list:

  A nested list with structure:
  `list(method = list(signature = data.frame(...)))`.

## Value

A data frame combining all the nested data frames, with added columns
`method` and `signature`.
