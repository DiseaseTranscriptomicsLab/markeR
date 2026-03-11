# Wrap Long Titles with Capital Letter Prioritization (when no symbols are nearby)

This function wraps long titles into multiple lines to fit a specified
width, prioritizing breaks at symbols like underscores, hyphens, and
colons when they are close to the wrap point. If no special symbol is
found nearby, the function will break the title at the first capital
letter. If neither is found, the title is broken at the specified width.

## Usage

``` r
wrap_title(title, width = 30)
```

## Arguments

- title:

  A character string representing the title to be wrapped.
  **(Required)**

- width:

  A numeric value specifying the maximum width for each line. The
  default is 30 characters. **(Optional)**

## Value

A character string with the title wrapped into multiple lines. Each line
will not exceed the specified width, with breaks prioritized by symbols
when nearby, and capital letters used only when no symbols are present.
