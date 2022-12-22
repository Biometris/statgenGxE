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
#' @name varCov
NULL

#' @rdname varCov
#' @keywords internal
createVarCov <- function(STA,
                         choice,
                         summary,
                         vcov,
                         criterion,
                         engine,
                         dat,
                         trait) {
  varCov <- structure(list(STA = STA,
                           choice = choice,
                           summary = summary,
                           vcov = vcov,
                           criterion = criterion,
                           engine = engine,
                           dat = dat,
                           trait = trait),
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
#' @param ... Not used.
#' @param title A character string used a title for the plot.
#' @param output Should the plot be output to the current device? If
#' \code{FALSE} only a ggplot object is invisibly returned.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("asreml", quietly = TRUE)) {
#'   ## Select the best variance-covariance model using asreml for modeling.
#'   geVarCov <- gxeVarCov(TD = TDMaize, trait = "yld", engine = "asreml")
#'
#'   ## Create a heatmap of the correlation matrix for the best model.
#'   plot(geVarCov)
#'   }
#' }
#'
#' @family varCov
#'
#' @export
plot.varCov <- function(x,
                        title = paste("Heatmap for model:", x$choice),
                        ...,
                        output = TRUE) {
  chkChar(title, len = 1, null = FALSE)
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
  p <- ggplot2::ggplot(data = meltedCorMatLow,
                       ggplot2::aes(x = .data[["trial1"]],
                                    y = .data[["trial2"]],
                                    fill = .data[["cor"]])) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_y_discrete(position = "right") +
    ## Create a gradient scale.
    ggplot2::scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                                  na.value = "grey", limit = c(-1, 1)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1,
                                                       size = 10, hjust = 1)) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ## Remove grid behind text output.
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::labs(title = title, x = "", y = "", fill = "correlation") +
    ## Fix coordinates to get a square sized plot.
    ggplot2::coord_fixed()
  if (output) {
    plot(p)
  }
  invisible(p)
}

#' Extract fitted values.
#'
#' Extract the fitted values for an object of class varCov.
#'
#' @param object An object of class varCov
#' @param ... Not used.
#'
#' @return A data.frame with fitted values.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("asreml", quietly = TRUE)) {
#'   ## Select the best variance-covariance model using asreml for modeling.
#'   geVarCov <- gxeVarCov(TD = TDMaize, trait = "yld", engine = "asreml")
#'   ## Extract fitted values from the model.
#'
#'   fitVarCov <- fitted(geVarCov)
#'   head(fitVarCov)
#'   }
#' }
#'
#' @family varCov
#'
#' @export
fitted.varCov <- function(object,
                          ...) {
  engine <- object$engine
  fitMod <- object$STA[[1]]$mFix[[1]]
  modDat <- object$dat
  if (engine == "lme4") {
    fittedValue <- predict(fitMod, newdata = modDat)
    fittedGeno <- cbind(modDat[c("trial", "genotype")], fittedValue)
  } else if (engine == "asreml") {
    fittedGeno <- predictAsreml(fitMod, classify = "trial:genotype",
                                TD = modDat)$pvals
    colnames(fittedGeno)[colnames(fittedGeno) == "predicted.value"] <-
      "fittedValue"
    colnames(fittedGeno)[colnames(fittedGeno) == "std.error"] <- "seFittedValue"
    fittedGeno <- fittedGeno[c("trial", "genotype", "fittedValue", "seFittedValue")]
  }
  return(fittedGeno)
}

#' Extract residuals.
#'
#' Extract the residuals for the best model.
#'
#' @param object An object of class varCov
#' @param ... Not used.
#'
#' @return A data.frame with residuals.
#'
#' @examples
#' \donttest{
#' ## Select the best variance-covariance model using asreml for modeling.
#' if (requireNamespace("asreml", quietly = TRUE)) {
#'   geVarCov <- gxeVarCov(TD = TDMaize, trait = "yld", engine = "asreml")
#'
#'   ## Extract residuals from the model.
#'   residVarCov <- residuals(geVarCov)
#'   head(residVarCov)
#'   }
#' }
#'
#' @family varCov
#'
#' @export
residuals.varCov <- function(object,
                             ...) {
  trait <- object$trait
  engine <- object$engine
  fittedGeno <- fitted(object)
  residGeno <- merge(fittedGeno, object$dat, by = c("trial", "genotype"))
  residGeno[["residual"]] <- residGeno[["fittedValue"]] - residGeno[[trait]]
  residGeno <- residGeno[c("trial", "genotype", "residual")]
  return(residGeno)
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
#' \donttest{
#' ## Select the best variance-covariance model using asreml for modeling.
#' if (requireNamespace("asreml", quietly = TRUE)) {
#'   geVarCov <- gxeVarCov(TD = TDMaize, trait = "yld", engine = "asreml")
#'
#'   ## Create a pdf report summarizing the results.
#'   report(geVarCov, outfile = tempfile(fileext = ".pdf"))
#'   }
#' }
#'
#' @family varCov
#'
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
