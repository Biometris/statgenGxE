#' S3 class varComp
#'
#' Function for creating objects of S3 class varComp.\cr
#' \code{\link{print}}, \code{\link{summary}}, \code{\link{plot}} and
#' \code{\link{report}} methods are available.
#'
#' @param SSA An object of class SSA, the best fitted model.
#' @param choice A character string indicating the best fitted model.
#' @param summary A data.frame with a summary of the fitted models.
#' @param vcov The covariance matrix of the best fitted model.
#' @param criterion A character string indicating the goodness-of-fit criterion
#' used for determinening the best model.
#' @param engine A character string containing the engine used for the analysis.
#'
#' @author Bart-Jan van Rossum
#'
#' @seealso \code{\link{plot.varComp}}, \code{\link{report.varComp}}
#'
#' @name varComp
NULL

#' @rdname varComp
#' @keywords internal
createVarComp <- function(SSA,
                          choice,
                          summary,
                          vcov,
                          criterion,
                          engine) {
  varComp <- structure(list(SSA = SSA,
                            choice = choice,
                            summary = summary,
                            vcov = vcov,
                            criterion = criterion,
                            engine = engine),
                       class = "varComp")
  attr(varComp, which = "timestamp") <- Sys.time()
  return(varComp)
}

#' @export
print.varComp <- function(x, ...) {
  cat(paste0("Best model: ", x$choice, ", based on ", x$criterion, ".\n\n"))

  x$summary
}

#' @export
summary.varComp <- function(object, ...) {
  print(object, ...)
}

#' Plot function for class varComp
#'
#' Function for plotting a heatmap of the correlation matrix for objects of
#' class varComp.
#'
#' @param x An object of class varComp
#' @param ... Not used
#' @param output Should the plot be output to the current device? If
#' \code{FALSE} only a list of ggplot objects is invisibly returned.
#'
#' @examples
#' \dontrun{
#' ## Select the best variance-covariance model using asreml for modeling.
#' geVarComp <- gxeVarComp(TD = TDMaize, trait = "yld", engine = "asreml")
#' ## Create a heatmap of the correlation matrix for the best model.
#' plot(geVarComp)
#' }
#'
#' @import stats
#' @export
plot.varComp <- function(x,
                         ...,
                         output = TRUE) {
  corMat <- cov2cor(x$vcov)
  PC1 <- princomp(corMat)$loadings[, 1]
  orderPC1 <- order(PC1)
  corMat <- corMat[orderPC1, orderPC1]
  varMat <- x$vcov[orderPC1, orderPC1]
  ## Melt variance and correlation matrices to get proper shape for ggplot.
  meltedCorMat <- reshape2::melt(corMat)
  meltedVarMat <- reshape2::melt(varMat)
  ## If trial names consist of only numbers melt converts them to numeric.
  ## This gives problems with plotting, so reconvert them to factor.
  if (is.numeric(meltedCorMat[["Var1"]])) {
    meltedCorMat[["Var1"]] <- factor(meltedCorMat[["Var1"]],
                                     levels = rownames(corMat))
    meltedCorMat[["Var2"]] <- factor(meltedCorMat[["Var2"]],
                                     levels = rownames(corMat))
  }
  if (is.numeric(meltedCorMat[["Var1"]])) {
    meltedVarMat[["Var1"]] <- factor(meltedVarMat[["Var1"]],
                                     levels = rownames(varMat))
    meltedVarMat[["Var2"]] <- factor(meltedVarMat[["Var2"]],
                                     levels = rownames(varMat))
  }
  ## Select bottom triangle for correlations and top for variances.
  meltedCorMatLow <- meltedCorMat[as.numeric(meltedCorMat$Var1) >
                                    as.numeric(meltedCorMat$Var2), ]
  meltedVarMatUp <- meltedVarMat[as.numeric(meltedVarMat$Var1) <=
                                   as.numeric(meltedVarMat$Var2), ]
  ## Round values for nicer display
  meltedVarMatUp$value <- round(meltedVarMatUp$value)
  p <- ggplot2::ggplot(data = meltedCorMatLow,
                  ggplot2::aes_string("Var1", "Var2", fill = "value")) +
    ggplot2::geom_tile(color = "white") +
    ## Discrete scales for x and y are needed to assure the diagonal is
    ## included in the plot. It is filled with variances later.
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ## Create a gradient scale.
    ggplot2::scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                                  na.value = "grey", limit = c(-1, 1)) +
    ## Var1 and Var2 have to be converted to numeric here to prevent ggplot
    ## from refactoring the data.
    ggplot2::geom_text(data = meltedVarMatUp,
                       ggplot2::aes_string("as.numeric(Var1)", "as.numeric(Var2)",
                         label = "value", size = "value")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1,
                                                       size = 10, hjust = 1)) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ## Remove grid behind text output.
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::ggtitle(paste("Heatmap for model:", x$choice)) +
    ggplot2::xlab("") + ggplot2::ylab("") +
    ggplot2::labs(fill = "correlation", size = "covariance") +
    ## Fix coordinates to get a square sized plot.
    ggplot2::coord_fixed()
  if (output) {
    plot(p)
  }
  invisible(p)
}

#' Report method for class varComp
#'
#' A pdf report will be created containing a summary of an object of class
#' varComp. Simultaneously the same report will be created as a tex
#' file.
#'
#' @inheritParams report.AMMI
#'
#' @param x An object of class varComp.
#'
#' @return A pdf and tex report.
#'
#' @examples
#' \dontrun{
#' ## Select the best variance-covariance model using asreml for modeling.
#' geVarComp <- gxeVarComp(TD = TDMaize, trait = "yld", engine = "asreml")
#' ## Create a pdf report summarizing the results.
#' report(geVarComp, outfile = "./testReports/reportVarComp.pdf")
#' }
#' @export
report.varComp <- function(x,
                           ...,
                           outfile = NULL) {
  createReport(x = x, reportName = "varCompReport.Rnw", outfile = outfile, ...)
}
