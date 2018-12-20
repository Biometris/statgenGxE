#' S3 class multiQTL
#'
#' Function for creating objects of S3 class multiQTL.\cr
#' \code{\link{summary}}, \code{\link{plot}} and \code{\link{report}}
#' methods are available.
#'
#' @param qtl A fitted multi QTL model.
#' @param QTLDet The object of class \code{\link{QTLDet}} used as base for
#' fitting the QTL model.
#' @param selection A character string indictating the type of selection used
#' for selecting the markers in the final model.
#' @param thr A numerical value indicating the threshold for dropping terms in
#' the backwards elemination process.
#'
#' @author Bart-Jan van Rossum
#'
#' @seealso \code{\link{plot.multiQTL}}, \code{\link{report.multiQTL}}
#'
#' @name multiQTL
NULL

#' @rdname multiQTL
#' @keywords internal
createMultiQTL <- function(qtl,
                           QTLDet,
                           selection,
                           thr) {
  multiQTL <- structure(list(qtl = qtl,
                             QTLDet = QTLDet,
                             selection = selection,
                             thr = thr),
                        class = "multiQTL")
  attr(multiQTL, which = "timestamp") <- Sys.time()
  return(multiQTL)
}

#' @export
print.multiQTL <- function(x, ...) {
  summary(x, ...)
}

#' @export
summary.multiQTL <- function(object, ...) {
  summary(object$qtl, ...)
}

#' Plot function for class multiQTL
#'
#' Plotting function for objects of class multiQTL. Plots the estimates of the
#' QTLs in the final model with their respective confidence intervals.
#'
#' @param x An object of class multiQTL
#' @param ... Further graphical parameters. Currently not used.
#' @param title A character string used as overall title for the plot.
#' @param output Should the plot be output to the current device? If
#' \code{FALSE} only a list of ggplot objects is invisibly returned.
#'
#' @examples
#' ## Read the data
#' F2 <- qtl::read.cross(format="csv",
#'                       file = system.file("extdata", "F2maize_geno.csv",
#'                                         package = "RAP"),
#'                       genotypes = c("AA", "AB", "BB"),
#'                       alleles = c("A", "B"), estimate.map = FALSE)
#' ## Perform QTL detection using simple interval mapping.
#' QTLDet <- QTLDetect(cross = F2, trait = "trait", type = "SIM")
#' ## Fit a multi QTL model.
#' multiFit <- multiQTLFit(QTLDet)
#' ## Plot the results.
#' plot(multiFit)
#'
#' @export
plot.multiQTL <- function(x,
                          ...,
                          title = "QTL estimates and confidence intervals",
                          output = TRUE) {
  ## Create plotData by extracting estimates and SE from summary.
  ## Remove intercept.
  summ <- summary(x$qtl)
  plotData <- as.data.frame(summ[["ests"]][-1, 1:2, drop = FALSE])
  ## Define as factor to prevent ggplot from plotting QTLs alphabetically.
  plotData$qtl <- factor(rownames(plotData), levels = rownames(plotData))
  p <- ggplot2::ggplot(plotData, ggplot2::aes_string(x = "qtl", y = "est")) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_errorbar(ggplot2::aes_string(ymin = "est - SE",
                                               ymax = "est + SE"),
                           width = 0.1) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                      vjust = 0.3, hjust = 1),
                   plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::ggtitle(title) +
    ggplot2::ylab("estimates") +
    ggplot2::geom_hline(yintercept = 0)
  if (output) {
    plot(p)
  }
  invisible(p)
}

#' Report method for class multiQTL
#'
#' A pdf report will be created containing a summary of a multiQTL analysis.
#' Simultaneously the same report will be created as a tex file.
#'
#' @inheritParams report.AMMI
#'
#' @param x An object of class multiQTL.
#'
#' @return A pdf and tex report.
#'
#' @examples
#' ## Read the data
#' F2 <- qtl::read.cross(format="csv",
#'                       file = system.file("extdata", "F2maize_geno.csv",
#'                                         package = "RAP"),
#'                       genotypes = c("AA", "AB", "BB"),
#'                       alleles = c("A", "B"), estimate.map = FALSE)
#' ## Perform QTL detection using simple interval mapping.
#' QTLDet <- QTLDetect(cross = F2, trait = "trait", type = "SIM")
#' ## Fit a multi QTL model.
#' multiFit <- multiQTLFit(QTLDet)
#' \dontrun{
#' ## Create a pdf report summarizing results.
#' report(multiFit, outfile = "./testReports/reportMultiQTLFit.pdf")
#' }
#'
#' @export
report.multiQTL <- function(x,
                            ...,
                            outfile = NULL) {
  createReport(x = x, reportName = "multiQTLReport.Rnw",
               outfile = outfile, ...)
}


