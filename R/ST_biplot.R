#' Draw a biplot on the current graphics device
#'
#' Draws a biplot on the current graphics device.
#' @param x A data frame or a matrix containing the data.
#' @param center  A logical value indicating whether the variables
#'  should be shifted to be zero centred. Alternately, a vector of
#'  length equal to the number of columns of \code{x} can be supplied.
#'  The value is passed to \code{scale}.
#' @param scale A logical value indicating whether the variables should
#'  be scaled to have unit variance before the analysis takes
#'  place. The default is \code{FALSE} for consistency with \code{\link{biplot.default}}, but
#'  in general scaling is advisable.  Alternatively, a vector of length
#'  equal to the number of columns of \code{x} can be supplied.  The
#'  value is passed to \code{\link{scale}}.
#' @param ... Additional arguments to be passed to \code{\link{biplot.default}}.
#'
#' @export

ST.biplot <- function(x,
                      center = TRUE,
                      scale = FALSE,
                      ...) {
  X <- as.matrix(x)
  gNames <- rownames(X)
  if (is.null(gNames)) {
    gNames <- 1:nrow(X)
  }
  nGNames <- sapply(X = as.character(gNames), FUN = nchar)
  if (mean(nGNames) <= 5) {
    rownames(X) <- gNames
  } else {
    rownames(X) <- 1:length(gNames)
  }
  pca <- prcomp(x = na.omit(X), center = center, scale. = scale)
  cump2 <- summary(pca)$importance[3, 2]
  propPC1 <- summary(pca)$importance[2, 1]
  propPC2 <- summary(pca)$importance[2, 2]
  biplot(pca, main = paste0("Principal components biplot (", round(cump2 * 100, 1), "%)"),
         xlab = paste0("PC1 (", round(propPC1 * 100, 1), "%)"),
         ylab = paste0("PC2 (", round(propPC2 * 100, 1), "%)") ,...)
}
