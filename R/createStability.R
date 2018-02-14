#' S3 class stability
#'
#' Function for creating objects of S3 class stability.\cr
#' \code{\link{print}}, \code{\link{summary}}, \code{\link{plot}} and
#' \code{\link{report}} methods are available.
#'
#' @param superiority a data.frame containing values for the cultivar-superiority
#' measure of Lin and Binns.
#' @param static a data.frame containing values for Shukla's stabilitye variance.
#' @param wricke a data.frame containing values for Wricke's ecovalence.
#' @param trait a character string indicating the trait that has been analyzed.
#'
#' @author Bart-Jan van Rossum
#'
#' @seealso \code{\link{plot.stability}}, \code{\link{report.stability}}
#'
#' @name stability
NULL

#' @rdname stability
#' @export
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
print.stability <- function(x, ...) {
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
summary.stability <- function(object, ...) {
  print(object, ...)
}

#' Plot Function for Class stability
#'
#' Function for creating plots of computed stability measures against the means.
#'
#' @param x an object of class stability.
#' @param ... other arguments to be passed on to underlying plot functions.
#'
#' @return Plots of stability measures against means.
#'
#' @import graphics grDevices
#' @export

plot.stability <- function(x,
                           ...) {
  nPlots <- sum(c(!is.null(x$superiority), !is.null(x$static), !is.null(x$wricke)))
  ## Prepare panels.
  if (nPlots > 1) {
    oldPar <- par(mfrow = c(nPlots - 1, 2), mar = c(4, 4, 2, 2),
                  oma = c(.5, .5, 2, .3))
    on.exit(par(oldPar))
  }
  ## Create superiority plot.
  if (!is.null(x$superiority)) {
    plot(x = x$superiority$mean, y = x$superiority$superiority,
         xlab = "Mean", ylab = "Cultivar superiority", ...)
  }
  ## Create static plot.
  if (!is.null(x$static)) {
    plot(x = x$static$mean, y = x$static$static,
         xlab = "Mean", ylab = "Static stability", ...)
  }
  ## Create Wricke plot.
  if (!is.null(x$wricke)) {
    plot(x = x$wricke$mean, y = x$wricke$wricke,
         xlab = "Mean", ylab = "Wricke's ecovalence", ...)
  }
  ## Add title.
  if (nPlots == 1) {
    title(paste('Stability coefficients for', x$trait))
  } else {
    mtext(paste('Stability coefficients for', x$trait), side = 3,
          outer = TRUE, cex = 1.3, font = 2)
  }
}

#' Report method for class stability
#'
#' A pdf report will be created containing a summary of stability model.
#' Simultaneously the same report will be created as a tex file.
#'
#' @param x an object of class stability.
#' @param ... further arguments passed on from other functions - not used yet.
#' @param outfile a character string, the name and location of the output .pdf and .tex
#' file for the report. If \code{NULL} a report will be created in the current working
#' directory.
#'
#' @export
report.stability <- function(x,
                             ...,
                             outfile = NULL) {
  createReport(x = x, reportName = "stabilityReport.Rnw",
               outfile = outfile, ...)
}



