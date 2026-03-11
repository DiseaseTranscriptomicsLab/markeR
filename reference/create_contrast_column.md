# Create Contrast Column in Metadata

This function extracts and processes contrast groups from a given
contrast string, then assigns contrast labels to a metadata subset based
on the variable of interest.

## Usage

``` r
create_contrast_column(metadata, variable_name, contrast)
```

## Arguments

- metadata:

  A data frame containing sample metadata.

- variable_name:

  A character string specifying the column name in `metadata` that
  represents the variable of interest.

- contrast:

  A character string representing the contrast in the form "(A + B) -
  (C + D)" (e.g.).

## Value

A subset of `metadata` with an added `cohentest` column, indicating
group membership based on the contrast.
