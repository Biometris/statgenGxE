#' Data Input ST
#'
#' Read a comma-separated values file obtained from a field trial study and create a data.frame from it.
#'
#' @param csvfile The name of the comma-separated values file which the data are to be read from.
#' @param factorNames A character vector of design factor names.
#' @param traitNames A character vector specifying trait names (and covariate names).
#' @param checkId (Optional) a string specifying the column name of the check ID(s).
#' @param env A string specifying the column name of the environment ID(s).
#' @param sep The field separator character. Values on each line of the file are separated
#' by this character. If sep = "" (the default for read.table) the separator is 'white space',
#' that is one or more spaces, tabs, newlines or carriage returns.
#' @param quote The set of quoting characters. To disable quoting altogether, use
#' quote = "". See \code{\link{scan}} for the behaviour on quotes embedded in quotes.
#' Quoting is only considered for columns read as character, which is all of them unless
#' colClasses is specified.
#' @param dec The character used in the file for decimal points.
#' @param fill Logical. If TRUE then in the case where the rows have unequal length, blank fields
#' are implicitly added. See details in \code{\link{read.table}}.
#' @param commentChar Character: a character vector of length one containing a single character
#' or an empty string. Use "" to turn off the interpretation of comments altogether.
#' @param rowSelect Character: a character vector of unique names to be selected from the
#' environment ID(s). If \code{NA}, all rows are included.
#' @param colSelect Character: a character vector of unique names to be selected for columns.
#' If \code{NA}, all columns are included.
#' @param naStrings A character vector of strings which are to be interpreted as \code{NA} values.
#' Blank fields are also considered to be missing values in logical, integer, numeric and
#' complex fields.
#' @param checkNames Logical. If \code{TRUE} then the names of the variables in the data.frame
#' are checked to ensure that they are syntactically valid variable names. If necessary
#' they are adjusted (by \code{\link{make.names}}) so that they are, and also to ensure that
#' there are no duplicates.
#' @param ... Further arguments to be passed to read.table.
#' @return A data.frame containing a representation of (a subset of) the data in the file.
#'
#' @seealso \code{\link{read.csv}} which this function wraps.
#'
#' @note a header is needed in the csv file which the data are to be read from.
#' @examples
#' myDat <- ST.read.csv(system.file("extdata", "SB_yield.csv", package = "RAP"),
#'                      factorNames = c("Env", "Genotype", "Rep", "Subblock", "Row", "Column"),
#'                      traitNames = "yield", env = "Env",rowSelect = c("HEAT06", "HEAT05"),
#'                      colSelect = c("Env", "Genotype", "Rep", "Row", "Column", "yield"))
#' str(myDat)
#'
#' @export
ST.read.csv <- function(csvfile,
                        factorNames,
                        traitNames,
                        checkId = NA,
                        env = NA,
                        sep = ",",
                        quote = "\"",
                        dec = ".",
                        fill = TRUE,
                        commentChar = "",
                        rowSelect = NA,
                        colSelect = NA,
                        naStrings = "*",
                        checkNames = FALSE,
                        ...) {
  # read the header lines as char strings
  hl <- readLines(csvfile, 1)
  # split headers up by 'sep'
  hl <- unlist(strsplit(hl, sep))
  hl <- gsub(pattern = quote, replacement = "", x = hl)
  nHl <- length(hl)
  if (nHl != length(unique(hl))) {
    stop("There are some duplicates of the names on the header.\n")
  }
  cls <- rep(NA, nHl)
  # the variables to be factorized
  if (is.na(checkId)) {
    fNames <- sapply(X = hl, FUN = function(x) {
      x %in% factorNames
    })
  } else {
    fNames <- sapply(X = hl, FUN = function(x) {
      x %in% c(factorNames, checkId)
    })
  }
  cls[fNames] <- "factor"
  # the variables to be numeric
  tNames <- sapply(X = hl, FUN = function(x) {
    x %in% traitNames
  })
  cls[tNames] <- "numeric"
  # columns to be selected...
  if (!all(is.na(colSelect))) {
    cls[!sapply(X = hl, FUN = function(x) {
      x %in% colSelect
    })] <- "NULL"
  }
  # use read.csv to read the data set
  myDat <- read.csv(csvfile, header = TRUE, sep = sep, quote = quote, dec = dec,
                    na.strings = naStrings, colClasses = cls, check.names=checkNames,
                    fill = fill, comment.char = commentChar, ...)
  # checks if the number of columns is correct
  if(all(is.na(colSelect))) {
    if (ncol(myDat) != nHl) {
      warning("CHECK: The length of column names on the header is not equal to the
              number of columns read in.\n")
    }
  } else {
    if (ncol(myDat) != length(colSelect)) {
      warning("CHECK: The number of selected columns is not equal to the number of columns
              read in.\n")
    }
  }
  # use data slicing to select rows according to environment id(s)
  if (env %in% hl) {
    if (!all(is.na(rowSelect))) {
      if (!all(is.na(colSelect))) {
        if (env %in% colSelect) {
          myDat <- myDat[myDat[[env]] %in% rowSelect, ]
          myDat <- droplevels(myDat)
        }
      } else {
        myDat <- myDat[myDat[[env]] %in% rowSelect,]
        myDat <- droplevels(myDat)
      }
      #for (dir in rowSelect) dir.create(dir)
      # check if the number of rows is correct
      if (length(unique(myDat[[env]])) != length(rowSelect)) {
        warning("CHECK: Either some selected row(s) are not read in\n or the selected
                row name(s) do not match the enviroment ID(s).\n")
      }
    }
  } else {
    # Add enviroment column
    message(paste0("Add an env column containing the repeated character string '", env, "'.\n"))
    env <- factor(rep(env, nrow(myDat)))
    myDat <- cbind(env, myDat)
  }
  return(myDat)
}
