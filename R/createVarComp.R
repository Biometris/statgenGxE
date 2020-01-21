#' S3 class varComp
#'
#' Function for creating objects of S3 class varComp.\cr
#' \code{\link{print}}, \code{\link{summary}}, \code{\link{plot}} and
#' \code{\link{report}} methods are available.
#'
#' @param STA An object of class STA, the best fitted model.
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
createVarComp <- function(STA,
                          choice,
                          summary,
                          vcov,
                          criterion,
                          engine) {
  varComp <- structure(list(STA = STA,
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
#' @export
plot.varComp <- function(x,
                         ...,
                         output = TRUE) {
  corMat <- cov2cor(x$vcov)
  PC1 <- princomp(corMat)$loadings[, 1]
  orderPC1 <- order(PC1)
  corMat <- corMat[orderPC1, orderPC1]
  ## Melt correlation matrix to get proper shape for ggplot.
  meltedCorMat <- reshape2::melt(corMat)
  ## If trial names consist of only numbers melt converts them to numeric.
  ## This gives problems with plotting, so reconvert them to factor.
  if (is.numeric(meltedCorMat[["Var1"]])) {
    meltedCorMat[["Var1"]] <- factor(meltedCorMat[["Var1"]],
                                     levels = rownames(corMat))
    meltedCorMat[["Var2"]] <- factor(meltedCorMat[["Var2"]],
                                     levels = rownames(corMat))
  }
  ## Select bottom triangle for correlations and top for variances.
  meltedCorMatLow <- meltedCorMat[as.numeric(meltedCorMat$Var1) >
                                    as.numeric(meltedCorMat$Var2), ]
  p <- ggplot(data = meltedCorMatLow,
              aes_string("Var1", "Var2", fill = "value")) +
    geom_tile(color = "white") +
    scale_y_discrete(position = "right") +
    ## Create a gradient scale.
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         na.value = "grey", limit = c(-1, 1)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                     size = 10, hjust = 1)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ## Remove grid behind text output.
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(title = paste("Heatmap for model:", x$choice), x = "", y = "",
         fill = "correlation") +
    ## Fix coordinates to get a square sized plot.
    coord_fixed()
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
  ## Checks.
  if (nchar(Sys.which("pdflatex")) == 0) {
    stop("An installation of LaTeX is required to create a pdf report.\n")
  }
  createReport(x = x, reportName = "varCompReport.Rnw", outfile = outfile, ...)
}
