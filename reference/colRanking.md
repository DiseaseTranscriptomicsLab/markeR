# Compute Independent Column-wise Ranks of Matrix Elements

This function computes the rank of each element in every column of a
numeric matrix independently. For each column, the smallest element
receives a rank of 1, the second smallest a rank of 2, and so on.

## Usage

``` r
colRanking(x, ties.method = "average")
```

## Arguments

- x:

  A numeric matrix.

- ties.method:

  A character string specifying the method used for tie-breaking.
  Options include `"average"`, `"first"`, `"random"`, `"max"`, or
  `"min"`. The default is `"average"`.

## Value

A numeric matrix of the same dimensions as `x` where each column
contains the ranks of the corresponding column's elements.
