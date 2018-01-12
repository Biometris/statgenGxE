#' S3 class varComp
#'
#' Function for creating objects of S3 class varComp.\cr
#' \code{\link{print}}, \code{\link{summary}}, \code{\link{plot}} and
#' \code{\link{report}} methods are available.
#'
#' @param model The best fitted model.
#' @param choice A character string indicating the best fitted model.
#' @param summary A data.frame with a summary of the fitted models.
#' @param vcov The covariance matrix of the best fitted model.
#' @param criterion A character string indicating the criterion used for
#' determinening the best model, either "AIC" or "BIC".
#' @param engine A character string containing the engine used for
#' the analysis.
#' @param x an \code{R} object
#'
#' @author Bart-Jan van Rossum
#'
#' @seealso \code{\link{plot.varComp}}, \code{\link{report.varComp}}
#'
#' @name varComp
NULL

#' @rdname varComp
#' @export
createVarComp <- function(model,
                          choice,
                          summary,
                          vcov,
                          criterion,
                          engine) {
  varComp <- structure(list(model = model,
                            choice = choice,
                            summary = summary,
                            vcov = vcov,
                            criterion = criterion,
                            engine = engine),
                       class = "varComp")
  attr(varComp, which = "timestamp") <- Sys.time()
  return(varComp)
}

#' @rdname varComp
#' @export
is.varComp <- function(x) {
  inherits(x, "varComp")
}

#' @export
print.varComp <- function(x, ...) {
  x$summary
}

#' @export
summary.varComp <- function(object, ...) {
  print(object, ...)
}

#' Plot Function for Class varComp
#'
#' Function for plotting a heatmap of the correlation matrix for objects of
#' class varComp.
#'
#' @param x an object of class varComp
#' @param ... not unused
#'
#' @import stats
#' @export
plot.varComp <- function(x, ...) {
  corMat <- cov2cor(x$vcov)
  meltedCorMat <- reshape2::melt(corMat)
  ggplot2::ggplot(data = meltedCorMat, ggplot2::aes_string("X1", "X2", fill = "value")) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, limit = c(-1, 1), space = "Lab") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1,
                                     size = 10, hjust = 1)) +
    ggplot2::xlab("") + ggplot2::ylab("") +
    ggplot2::coord_fixed()
}

#' Report method for class varComp
#'
#' A pdf report will be created containing a summary of selection of the
#' best variance model.
#' Simultaneously the same report will be created as a tex file.
#'
#' @param x an object of class varComp.
#' @param ... further arguments passed on from other functions - not used yet.
#' @param outfile a character string, the name and location of the output .pdf and .tex
#' file for the report. If \code{NULL} a report will be created in the current working
#' directory.
#'
#' @export
report.varComp <- function(x,
                           ...,
                           outfile = NULL) {
  createReport(x = x, reportName = "varCompReport.Rnw",
               outfile = outfile, ...)
}



