#' Wrap Long Titles with Capital Letter Prioritization (when no symbols are nearby)
#'
#' This function wraps long titles into multiple lines to fit a specified width, prioritizing breaks at symbols like
#' underscores, hyphens, and colons when they are close to the wrap point. If no special symbol is found nearby, the function
#' will break the title at the first capital letter. If neither is found, the title is broken at the specified width.
#'
#' @param title A character string representing the title to be wrapped. **(Required)**
#' @param width A numeric value specifying the maximum width for each line. The default is 30 characters. **(Optional)**
#'
#' @return A character string with the title wrapped into multiple lines. Each line will not exceed the specified width,
#'   with breaks prioritized by symbols when nearby, and capital letters used only when no symbols are present.
#'
#' @keywords internal
wrap_title <- function(title, width = 30) {
  if (nchar(title) <= width) {
    return(title)  # No need to wrap if it fits
  }

  wrapped_title <- ""
  while (nchar(title) > width) {
    # Find positions of capital letters and symbols near the wrap point
    capital_pos <- gregexpr("[A-Z]", title)[[1]]
    symbol_pos <- gregexpr("(_|-|:|\\+|\\\\|/|\\*|\\.|,|;|\\?|!)", title)[[1]]

    # Check for symbol breaks within the last few characters (width - 5 to width)
    valid_symbol_breaks <- symbol_pos[symbol_pos >= (width - 5) & symbol_pos <= width]

    if (length(valid_symbol_breaks) > 0) {
      # If a suitable symbol is found, break at the first valid symbol
      break_at <- valid_symbol_breaks[1]
    } else {
      # If no suitable symbol, look for capital letters within the same range
      valid_capital_breaks <- capital_pos[capital_pos >= (width - 5) & capital_pos <= width]

      if (length(valid_capital_breaks) > 0) {
        # If a capital letter is found, break just before the capital letter
        break_at <- valid_capital_breaks[1] - 1
      } else {
        # If no suitable symbol or capital letter, break at width
        break_at <- width
      }
    }

    # Append the wrapped line
    wrapped_title <- paste0(wrapped_title, substr(title, 1, break_at), "\n")

    # Update title with the remaining text after the break
    title <- substr(title, break_at + 1, nchar(title))
  }

  # Add the remaining part of the title
  wrapped_title <- paste0(wrapped_title, title)

  return(wrapped_title)
}
