.onAttach <- function(libname, pkgname) {
  gg_version <- tryCatch(
    utils::packageVersion("ggplot2"),
    error = function(e) NULL
  )
  
  if (!is.null(gg_version) && gg_version > "3.5.2") {
    packageStartupMessage(
      "Warning: markeR has been tested with ggplot2 <= 3.5.2. ",
      "Using newer versions may cause incompatibilities."
    )
  }
}
