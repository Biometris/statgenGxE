#' @keywords internal
"_PACKAGE"
#'
#' @import stats
#' @importFrom graphics plot
#' @importFrom xtable xtable
#' @importFrom utils hasName
#' @importFrom statgenSTA createTD
NULL

statgenDefaultOptions <- list(
  ##palette.colors(palette = "Dark 2")
  statgen.genoColors = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
                         "#E6AB02", "#A6761D", "#666666"),
  ##palette.colors(palette = "Alphabet")
  statgen.trialColors = c("#AA0DFE", "#3283FE", "#85660D", "#782AB6", "#565656",
                          "#1C8356", "#16FF32", "#F7E1A0", "#E2E2E2", "#1CBE4F",
                          "#C4451C", "#DEA0FD", "#FE00FA", "#325A9B", "#FEAF16",
                          "#F8A19F", "#90AD1C", "#F6222E", "#1CFFCE", "#2ED9FF",
                          "#B10DA1", "#C075A6", "#FC1CBF", "#B00068", "#FBE426",
                          "#FA0087")
)

.onLoad <- function(libname, pkgname) {
  op <- options()
  toset <- !(names(statgenDefaultOptions) %in% names(op))
  if (any(toset)) {
    options(statgenDefaultOptions[toset])
  }

  return(invisible(NULL))
}
