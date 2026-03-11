# Compute Cohen\\s d Effect Size

Computes the absolute Cohen\\s d effect size between two numeric
vectors. This function returns the absolute value of the difference in
means divided by the pooled standard deviation.

## Usage

``` r
cohen_d(x, y)
```

## Arguments

- x:

  A numeric vector representing the values for group 1.

- y:

  A numeric vector representing the values for group 2.

## Value

A numeric value representing Cohen\\s d. Returns NA if either group has
fewer than two observations or if the pooled standard deviation is zero.
