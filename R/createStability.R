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
#' @name stability
NULL

#' @rdname stability
#' @keywords internal
createStability <- function(superiority = NULL,
                            static = NULL,
                            wricke = NULL,
                            TD,
                            trait) {
  stability <- structure(list(superiority = superiority,
                              static = static,
                              wricke = wricke,
                              TD = TD,
                              trait = trait),
                         class = "stability")
  attr(stability, which = "timestamp") <- Sys.time()
  return(stability)
}

#' @export
print.stability <- function(x,
                            ...,
                            pctGeno = 10) {
  if (!is.null(x$superiority)) {
    cat("\nCultivar-superiority measure (Top", pctGeno, "% genotypes)\n")
    print(x$superiority[1:(ceiling(nrow(x$superiority) / 100 * pctGeno)), ],
          row.names = FALSE)
  }
  if (!is.null(x$static)) {
    cat("\nStatic stability (Top", pctGeno, "% genotypes)\n")
    print(x$static[1:(ceiling(nrow(x$static) / 100 * pctGeno)), ],
          row.names = FALSE)
  }
  if (!is.null(x$wricke)) {
    cat("\nWricke's ecovalence (Top", pctGeno, "% genotypes)\n")
    print(x$wricke[1:(ceiling(nrow(x$wricke) / 100 * pctGeno)), ],
          row.names = FALSE)
  }
}

#' @export
summary.stability <- function(object,
                              ...) {
  print(object, ...)
}

#' Plot function for class stability
#'
#' Function for creating scatter plots of the square roots of the computed
#' stability measures against the means.
#'
#' @param x An object of class stability.
#' @param ... Not used.
#' @param colorGenoBy A character string indicating a column in the \code{TD}
#' used as input for the stability analysis by which the genotypes should be
#' colored. If \code{NULL} all genotypes will be colored black.
#' @param colGeno A character vector with plot colors for the genotypes. A
#' single color when \code{colorGenoBy = NULL}, a vector of colors otherwise.
#' @param title A character string used a title for the plot.
#' @param output Should the plot be output to the current device? If
#' \code{FALSE} only a list of ggplot objects is invisibly returned.
#'
#' @return A list of ggplot object is invisibly returned.
#'
#' @examples
#' ## Compute three stability measures for TDMaize.
#' geStab <- gxeStability(TD = TDMaize, trait = "yld")
#'
#' ## Create scatter plots of the computed stability measures against the means.
#' plot(geStab)
#'
#' @family stability
#'
#' @importFrom grDevices topo.colors
#' @export
plot.stability <- function(x,
                           ...,
                           colorGenoBy = NULL,
                           colGeno = NULL,
                           title = paste("Stability coefficients for", x$trait),
                           output = TRUE) {
  chkChar(title, len = 1, null = FALSE)
  TDTot <- do.call(rbind, x$TD)
  if (!is.null(colorGenoBy)) {
    chkCol(colorGenoBy, TDTot)
  }
  chkChar(colGeno)
  genoDat <- unique(TDTot[c("genotype", colorGenoBy)])
  if (!is.null(colorGenoBy)) {
    if (!is.factor(genoDat[[colorGenoBy]])) {
      genoDat[[colorGenoBy]] <- as.factor(genoDat[[colorGenoBy]])
    }
    genoDat <- genoDat[order(genoDat[[colorGenoBy]]), ]
  } else {
    genoDat[[".colorGenoBy"]] <- factor(1)
    colorGenoBy <- ".colorGenoBy"
  }
  nColGeno <- nlevels(genoDat[[colorGenoBy]])
  if (length(colGeno) == 0) {
    ## Defaults to black for one color for genotypes.
    ## For more than one colors from statgen.genoColors are used.
    ## Fall back to topo.colors if number of colors in option is too small.
    if (nColGeno == 1) {
      colGeno <- "black"
    } else if (length(getOption("statgen.genoColors")) >= nColGeno) {
      colGeno <- getOption("statgen.genoColors")[1:nColGeno]
    } else {
      colGeno <- topo.colors(nColGeno)
    }
  } else {
    nColGenoArg <- length(colGeno)
    if (nColGenoArg != nColGeno) {
      stop("Number of colors provided doesn't match number of genotype groups:",
           "\n", nColGenoArg, " colors provided, ", nColGeno,
           " groups in data.\n")
    }
  }
  plots <- vector(mode = "list")
  if (!is.null(x$superiority)) {
    ## Create superiority plot.
    supDat <- merge(x$superiority, genoDat, by.x = "Genotype", by.y = "genotype")
    plots$p1 <- ggplot2::ggplot(data = supDat,
                                ggplot2::aes(x = .data[["Mean"]],
                                             y = sqrt(.data[["Superiority"]]),
                                             color = .data[[colorGenoBy]])) +
      ggplot2::geom_point() +
      ggplot2::scale_color_manual(values = colGeno) +
      ggplot2::labs(x = "Mean", y = "Square root of\n Cultivar superiority")
  }
  if (!is.null(x$static)) {
    ## Create static plot.
    statDat <- merge(x$static, genoDat, by.x = "Genotype", by.y = "genotype")
    plots$p2 <- ggplot2::ggplot(data = statDat,
                                ggplot2::aes(x = "Mean",
                                             y = sqrt(.data[["Static"]]),
                                             color = .data[[colorGenoBy]])) +
      ggplot2::geom_point() +
      ggplot2::scale_color_manual(values = colGeno) +
      ggplot2::labs(x = "Mean", y = "Square root of\n Static stability")
  }
  if (!is.null(x$wricke)) {
    ## Create Wricke plot.
    wrickeDat <- merge(x$wricke, genoDat, by.x = "Genotype", by.y = "genotype")
    plots$p3 <- ggplot2::ggplot(data = wrickeDat,
                                ggplot2::aes(x = "Mean",
                                             y = sqrt(.data[["Wricke"]]),
                                             color = .data[[colorGenoBy]])) +
      ggplot2::geom_point() +
      ggplot2::scale_color_manual(values = colGeno) +
      ggplot2::labs(x = "Mean", y = "Square root of\n Wricke's ecovalence")
  }
  ## Construct legend.
  if (colorGenoBy != ".colorGenoBy") {
    ## Build plot to extract legend.
    ## Legend is always the same for all plots, take it from the first plot.
    p1Gtable <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(plots[[1]]))
    legendPos <- sapply(X = p1Gtable$grobs, FUN = `[[`, "name") == "guide-box"
    legend <- p1Gtable$grobs[[which(legendPos)]]
  } else {
    legend <- NULL
  }
  ## Now remove the legend from all the plots.
  for (i in seq_along(plots)) {
    plots[[i]] <- plots[[i]] + ggplot2::theme(legend.position = "none")
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

  if (output) {
    ## grid.arrange automatically plots the results.
    tot <- gridExtra::grid.arrange(tot, right = legend, top = title)
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
#' \donttest{
#' ## Create a .pdf report summarizing the stability measures.
#' report(geStab, outfile = tempfile(fileext = ".pdf"))
#' }
#'
#' @family stability
#'
#' @export
report.stability <- function(x,
                             ...,
                             outfile = NULL) {
  ## Checks.
  if (nchar(Sys.which("pdflatex")) == 0) {
    stop("An installation of LaTeX is required to create a pdf report.\n")
  }
  createReport(x = x, reportName = "stabilityReport.Rnw",
               reportPackage = "statgenGxE", outfile = outfile, ...)
}
