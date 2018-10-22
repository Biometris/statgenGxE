#' S3 class QTLDet
#'
#' Function for creating objects of S3 class QTLDet.\cr
#' \code{\link{print}}, \code{\link{summary}}, \code{\link{plot}} and
#' \code{\link{report}} methods are available.
#'
#' @param scores A data.frame containing the lod scores.
#' @param peaks A data.frame containing the peaks found.
#' @param type A character string indicating the type of QTL detection performed.
#' @param cross An object of class cross in the \code{qtl} package.
#' @param trait A character string indicating the trait for which the analysis
#' is done.
#' @param info A list containing information on the settings used for
#' QTL detection, i.e. step, threshold and window.
#'
#' @author Bart-Jan van Rossum
#'
#' @seealso \code{\link{plot.QTLDet}}, \code{\link{report.QTLDet}}
#'
#' @name QTLDet
NULL

#' @rdname QTLDet
#' @export
createQTLDet <- function(scores,
                         peaks,
                         type,
                         cross,
                         trait,
                         info) {
  QTLDet <- structure(list(scores = scores,
                           peaks = peaks,
                           type = type,
                           cross = cross,
                           trait = trait,
                           info = info),
                      class = "QTLDet")
  attr(QTLDet, which = "timestamp") <- Sys.time()
  return(QTLDet)
}

#' @export
print.QTLDet <- function(x, ...) {
  if (nrow(x$peaks) != 0) {
    cat("Peaks\n")
    cat("=====\n")
    print(x$peaks)
  } else {
    cat("No peaks detected")
  }
}

#' @export
summary.QTLDet <- function(object,
                           ...) {
  print(object, ...)
}

#' Plot function for class QTLDet
#'
#' Function for creating a manhattan plot for objects of class QTLDet.
#'
#' @param x An object of class QTLDet
#' @param ... Not used
#' @param yLim A numerical value for the upper limit of the y-axis. If
#' \code{NA} the limit is determined based on the data.
#' @param title A character value for the title of the plot. If not supplied
#' the trait used in the qTL Detection is used as plot title.
#' @param output Should the plot be output to the current device? If
#' \code{FALSE} only a list of ggplot objects is invisibly returned.
#'
#' @import graphics grDevices
#' @keywords internal
plot.QTLDet <- function(x,
                        ...,
                        yLim = NA,
                        title = x$trait,
                        output = TRUE) {
  if (!is.na(yLim) && (!is.numeric(yLim) || length(yLim) > 1 || yLim < 0)) {
    stop("yLim should be NA or a single positive numerical value.\n")
  }
  if (!is.character(title) || length(title) > 1) {
    stop("title should be a single character string.\n")
  }
  p <- ggplot2::ggplot(data = x$scores,
                       ggplot2::aes_string(x = "pos", y = "lod")) +
    ggplot2::geom_line() +
    ## Add grid for plotting per chromosome. Scales and space free for x to
    ## get proper width for all chromosomes.
    ggplot2::facet_grid(rows = formula("~chr"), scales = "free_x",
                        space = "free_x") +
    ## Set continuous scales to start axes at (0,0).
    ggplot2::scale_y_continuous(limits = c(0, yLim), expand = c(0, 0)) +
    ggplot2::scale_x_continuous(limits = c(0, NA), expand = c(0, 0)) +
    ## Add cartesian coordinate system to get access to clip off option.
    ## Needed for adding peaks.
    ggplot2::coord_cartesian(clip = "off") +
    ## Set axis en plot titles.
    ggplot2::labs(x = "Chromosome", y = "LOD") +
    ggplot2::ggtitle(label = title) +
    ggplot2::theme(plot.title =  ggplot2::element_text(hjust = 0.5))
  if (x$type == "CIM") {
    if (nrow(x$peaks) > 0) {
      ## Add peaks as red dots on x-axis.
      p <- p + ggplot2::geom_point(data = x$peaks,
                                   ggplot2::aes_string(y = 0), color = "red")
    }
  }
  if (output) {
    plot(p)
  }
  invisible(p)
}

#' Report method for class QTLDet
#'
#' A pdf report will be created containing a summary of a QTLDet analysis.
#' Simultaneously the same report will be created as a tex file.
#'
#' @inheritParams report.AMMI
#'
#' @param x An object of class QTLDet.
#'
#' @return A pdf and tex report.
#'
#' @examples
#' F2 <- qtl::read.cross(format="csv",
#'                       file = system.file("extdata",
#'                                         "F2_maize_practical3_ex2.csv",
#'                                         package = "RAP"),
#'                       genotypes = c("AA", "AB", "BB"),
#'                       alleles = c("A", "B"), estimate.map = FALSE)
#' ## Perform a composite interval mapping for detecting QTLs.
#' QTLDet <- QTLDetect(cross = F2, trait = "trait", type = "CIM")
#' \dontrun{
#' ## Create a pdf report summarizing the results.
#' report(QTLDet, outfile = "./testReports/reportQTLDectection.pdf")
#' }
#'
#' @export
report.QTLDet <- function(x,
                          ...,
                          outfile = NULL) {
  createReport(x = x, reportName = "QTLDetReport.Rnw",
               outfile = outfile, ...)
}


