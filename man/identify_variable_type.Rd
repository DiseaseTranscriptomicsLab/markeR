% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VariableAssociation.R
\name{identify_variable_type}
\alias{identify_variable_type}
\title{Identify Variable Types}
\usage{
identify_variable_type(df, cols = NULL)
}
\arguments{
\item{df}{A data frame containing the variables to classify.}

\item{cols}{A character vector of column names to consider.}
}
\value{
A named character vector where names correspond to column names
and values indicate the variable type: "Numeric", "Categorical Bin", or "Categorical Multi".
}
\description{
Determines the type of each variable in a given data frame.
Variables are classified as "Numeric", "Categorical Bin" (binary categorical),
or "Categorical Multi" (multi-level categorical). Warnings are issued if
categorical variables have more than 10 unique values.
}
\examples{
\dontrun{
df <- data.frame(
  age = c(25, 30, 35, 40),
  gender = c("Male", "Female", "Female", "Male"),
  score = c(80, 85, 90, 95)
)
identify_variable_type(df)
}
}
\keyword{internal}
