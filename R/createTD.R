#' S3 class TD
#'
#' Function for creating objects of S3 class TD (Trial Data).
#'
#' \code{\link{print}} and \code{\link{summary}} methods are available.
#'
#' @param data a data.frame containing trial data with a least a column for
#' genotype and environment.
#' @param genotype a character string indicating the column in \code{data} that
#' contains genotypes.
#' @param env a character string indicating the column in \code{data} that
#' contains environments.
#' @param megaEnv an optional character string indicating the column in \code{data} that
#' contains megaEnvironments as constructed by \code{\link{GE.megaEnvironment}}.
#' @param year an optional character string indicating the column in \code{data} that
#' contains years.
#' @param rep an optional character string indicating the column in \code{data} that
#' contains replicates.
#' @param subBlock an optional character string indicating the column in \code{data} that
#' contains sub blocks.
#' @param row an optional character string indicating the column in \code{data} that
#' contains field rows.
#' @param col an optional character string indicating the column in \code{data} that
#' contains field columns.

#' @param x an \code{R} object
#'
#' @author Bart-Jan van Rossum
#'
#' @seealso \code{\link{summary.TD}}
#'
#' @name TD
NULL

#' @rdname TD
#' @export
createTD <- function(data,
                     genotype = "genotype",
                     env = "env",
                     megaEnv = NULL,
                     year = NULL,
                     rep = NULL,
                     subBlock = NULL,
                     row = NULL,
                     col = NULL) {
  cols <- colnames(data)
  cols[cols == genotype] <- "genotype"
  cols[cols == env] <- "env"
  cols[cols == megaEnv] <- "megaEnv"
  cols[cols == year] <- "year"
  cols[cols == rep] <- "rep"
  cols[cols == subBlock] <- "subBlock"
  cols[cols == row] <- "row"
  cols[cols == col] <- "col"
  colnames(data) <- cols
  TD <- structure(data,
                  class = c("TD", "data.frame"))
  return(TD)
}

#' @rdname TD
#' @export
is.TD <- function(x) {
  inherits(x, "TD")
}




