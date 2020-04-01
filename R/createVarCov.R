#' S3 class varCov
#'
#' Function for creating objects of S3 class varCov.\cr
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
#' @seealso \code{\link{plot.varCov}}, \code{\link{report.varCov}}
#'
#' @name varCov
NULL

#' @rdname varCov
#' @keywords internal
createVarCov <- function(STA,
                          choice,
                          summary,
                          vcov,
                          criterion,
                          engine) {
  varCov <- structure(list(STA = STA,
                            choice = choice,
                            summary = summary,
                            vcov = vcov,
                            criterion = criterion,
                            engine = engine),
                       class = "varCov")
  attr(varCov, which = "timestamp") <- Sys.time()
  return(varCov)
}

#' @export
print.varCov <- function(x, ...) {
  cat(paste0("Best model: ", x$choice, ", based on ", x$criterion, ".\n\n"))

  x$summary
}

#' @export
summary.varCov <- function(object, ...) {
  print(object, ...)
}

#' Plot function for class varCov
#'
#' Function for plotting a heatmap of the correlation matrix for objects of
#' class varCov.
#'
#' @param x An object of class varCov
#' @param ... Not used
#' @param output Should the plot be output to the current device? If
#' \code{FALSE} only a list of ggplot objects is invisibly returned.
#'
#' @examples
#' \dontrun{
#' ## Select the best variance-covariance model using asreml for modeling.
#' geVarCov <- gxeVarCov(TD = TDMaize, trait = "yld", engine = "asreml")
#' ## Create a heatmap of the correlation matrix for the best model.
#' plot(geVarCov)
#' }
#'
#' @export
plot.varCov <- function(x,
                         ...,
                         output = TRUE) {
  corMat <- cov2cor(x$vcov)
  PC1 <- princomp(corMat)$loadings[, 1]
  orderPC1 <- order(PC1)
  corMat <- corMat[orderPC1, orderPC1]
  ## Convert corMat to data.frame to prevent crash when reshaping.
  corMat <- as.data.frame(corMat)
  ## Convert correlation matrix to long format for ggplot.
  meltedCorMat <- reshape(corMat, direction = "long",
                          varying = list(genotype = colnames(corMat)),
                          ids = rownames(corMat), idvar = "trial1",
                          times = colnames(corMat), timevar = "trial2",
                          v.names = "cor")
  ## Reshape converts trial columns to character.
  ## This gives problems with plotting, so reconvert them to factor.
  meltedCorMat[["trial1"]] <- factor(meltedCorMat[["trial1"]],
                                     levels = rownames(corMat))
  meltedCorMat[["trial2"]] <- factor(meltedCorMat[["trial2"]],
                                     levels = rownames(corMat))
  ## Select bottom right triangle for correlations and top for variances.
  meltedCorMatLow <- meltedCorMat[as.numeric(meltedCorMat[["trial1"]]) >
                                    as.numeric(meltedCorMat[["trial2"]]), ]
  p <- ggplot(data = meltedCorMatLow,
              aes_string("trial1", "trial2", fill = "cor")) +
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

#' Report method for class varCov
#'
#' A pdf report will be created containing a summary of an object of class
#' varCov. Simultaneously the same report will be created as a tex
#' file.
#'
#' @inheritParams report.AMMI
#'
#' @param x An object of class varCov.
#'
#' @return A pdf and tex report.
#'
#' @examples
#' \dontrun{
#' ## Select the best variance-covariance model using asreml for modeling.
#' geVarCov <- gxeVarCov(TD = TDMaize, trait = "yld", engine = "asreml")
#' ## Create a pdf report summarizing the results.
#' report(geVarCov, outfile = "./testReports/reportVarCov.pdf")
#' }
#' @export
report.varCov <- function(x,
                           ...,
                           outfile = NULL) {
  ## Checks.
  if (nchar(Sys.which("pdflatex")) == 0) {
    stop("An installation of LaTeX is required to create a pdf report.\n")
  }
  createReport(x = x, reportName = "varCovReport.Rnw", outfile = outfile,
               reportPackage = "statgenGxE", ...)
}