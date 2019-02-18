#' S3 class stability
#'
#' Function for creating objects of S3 class stability.\cr
#' \code{\link{print}}, \code{\link{summary}}, \code{\link{plot}} and
#' \code{\link{report}} methods are available.
#'
#' @param superiority A data.frame containing values for the
#' cultivar-superiority measure of Lin and Binns.
#' @param static A data.frame containing values for Shukla's stability variance.
#' @param wricke A data.frame containing values for Wricke's ecovalence.
#' @param trait A character string indicating the trait that has been analyzed.
#'
#' @author Bart-Jan van Rossum
#'
#' @seealso \code{\link{plot.stability}}, \code{\link{report.stability}}
#'
#' @name stability
NULL

#' @rdname stability
#' @keywords internal
createstability <- function(superiority = NULL,
                            static = NULL,
                            wricke = NULL,
                            trait) {
  stability <- structure(list(superiority = superiority,
                              static = static,
                              wricke = wricke,
                              trait = trait),
                         class = "stability")
  attr(stability, which = "timestamp") <- Sys.time()
  return(stability)
}

#' @export
print.stability <- function(x,
                            ...) {
  if (!is.null(x$superiority)) {
    cat("\nCultivar-superiority measure (Top 10% genotypes)\n")
    print(x$superiority[1:(ceiling(nrow(x$superiority) / 10)), ],
          row.names = FALSE)
  }
  if (!is.null(x$static)) {
    cat("\nStatic stability (Top 10% genotypes)\n")
    print(x$static[1:(ceiling(nrow(x$static) / 10)), ], row.names = FALSE)
  }
  if (!is.null(x$wricke)) {
    cat("\nWricke's ecovalence (Top 10% genotypes)\n")
    print(x$wricke[1:(ceiling(nrow(x$wricke) / 10)), ], row.names = FALSE)
  }
}

#' @export
summary.stability <- function(object,
                              ...) {
  print(object, ...)
}

#' Plot function for class stability
#'
#' Function for creating scatter plots of computed stability measures against
#' the means.
#'
#' @param x An object of class stability.
#' @param ... Further arguments to be passed on to underlying plot functions.
#' @param output Should the plot be output to the current device? If
#' \code{FALSE} only a list of ggplot objects is invisibly returned.
#'
#' @return Plots of stability measures against means.
#'
#' @examples
#' ## Compute three stability measures for TDMaize.
#' geStab <- gxeStability(TD = TDMaize, trait = "yld")
#' ## Create scatter plots of the computed stability measures against the means.
#' plot(geStab)
#'
#' @import graphics grDevices
#' @export
plot.stability <- function(x,
                           ...,
                           output = TRUE) {
  dotArgs <- list(...)
  ## Add and overwrite args with custom args from ...
  fixedArgs <- c("x", "y", "xlab", "ylab", "main")
  plots <- vector(mode = "list")
  if (!is.null(x$superiority)) {
    ## Create superiority plot.
    plots$p1 <- ggplot2::ggplot(data = x$superiority,
                                ggplot2::aes_string(x = "mean",
                                                    y = "superiority")) +
      ggplot2::geom_point(col = "blue", shape = 1) +
      ggplot2::labs(x = "Mean", y = "Cultivar superiority")
  }
  if (!is.null(x$static)) {
    ## Create static plot.
    plots$p2 <- ggplot2::ggplot(data = x$static,
                                ggplot2::aes_string(x = "mean",
                                                    y = "static")) +
      ggplot2::geom_point(col = "blue", shape = 1) +
      ggplot2::labs(x = "Mean", y = "Static stability")
  }
  if (!is.null(x$wricke)) {
    ## Create Wricke plot.
    plots$p3 <- ggplot2::ggplot(data = x$wricke,
                                ggplot2::aes_string(x = "mean",
                                                    y = "wricke")) +
      ggplot2::geom_point(col = "blue", shape = 1) +
      ggplot2::labs(x = "Mean", y = "Wricke's ecovalence")
  }
  if (length(plots) == 3) {
    ## Create empty plot for bottom right grid position.
    plots$p4 <- ggplot2::ggplot() +
      ggplot2::theme(panel.background = ggplot2::element_blank())
  }
  ## Convert plots to grob for outlining of axes.
  plotsGr <- lapply(X = plots, FUN = ggplot2::ggplotGrob)
  if (length(plotsGr) > 1) {
    ## At least two plots -> 1 row.
    tot <- gridExtra::gtable_cbind(plotsGr[[1]], plotsGr[[2]])
  } else {
    ## Only 1 plot to be made.
    tot <- plotsGr[[1]]
  }
  if (length(plotsGr) > 2) {
    ## 3 Plots, so empty plot has been created to make 4 plots.
    ## First create second row then add to first row.
    r2 <- gridExtra::gtable_cbind(plotsGr[[3]], plotsGr[[4]])
    tot <- gridExtra::gtable_rbind(tot, r2)
  }
  ## Construct title.
  if (!is.null(dotArgs$title)) {
    title <- dotArgs$title
  } else {
    title <- paste("Stability coefficients for", x$trait)
  }
  if (output) {
    ## grid.arrange automatically plots the results.
    tot <- gridExtra::grid.arrange(tot, top = title)
  }
  invisible(plots)
}

#' Report method for class stability
#'
#' A pdf report will be created containing a summary of an object of class
#' stability. Simultaneously the same report will be created as a tex file.
#'
#' @inheritParams report.AMMI
#'
#' @param x An object of class stability.
#'
#' @return A pdf and tex report.
#'
#' @examples
#' ## Compute three stability measures for TDMaize.
#' geStab <- gxeStability(TD = TDMaize, trait = "yld")
#' \dontrun{
#' ## Create a .pdf report summarizing the stability measures.
#' report(geStab, outfile = "./testReports/reportStability.pdf")
#' }
#'
#' @export
report.stability <- function(x,
                             ...,
                             outfile = NULL) {
  createReport(x = x, reportName = "stabilityReport.Rnw",
               outfile = outfile, ...)
}
