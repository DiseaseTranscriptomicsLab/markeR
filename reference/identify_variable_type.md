# Identify Variable Types

Determines the type of each variable in a given data frame. Variables
are classified as "Numeric", "Categorical Bin" (binary categorical), or
"Categorical Multi" (multi-level categorical). Warnings are issued if
categorical variables have more than 10 unique values.

## Usage

``` r
identify_variable_type(df, cols = NULL)
```

## Arguments

- df:

  A data frame containing the variables to classify.

- cols:

  A character vector of column names to consider.

## Value

A named character vector where names correspond to column names and
values indicate the variable type: "Numeric", "Categorical Bin", or
"Categorical Multi".

## Examples

``` r
if (FALSE) { # \dontrun{
df <- data.frame(
  age = c(25, 30, 35, 40),
  gender = c("Male", "Female", "Female", "Male"),
  score = c(80, 85, 90, 95)
)
identify_variable_type(df)
} # }
```
