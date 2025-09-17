.onAttach <- function(libname, pkgname) {
  if (packageVersion("ggplot2") > "3.5.2") {
    packageStartupMessage(
      "Warning: markeR has been tested with ggplot2 <= 3.5.2. ",
      "Using newer versions may cause incompatibilities."
    )
  }
}
