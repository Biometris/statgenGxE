#' Generates an HTML table
#'
#' Write a matrix or data.frame as an HTML table.
#'
#' @param x A matrix or data.frame object.
#' @param rowNames A logical specifying if the row names are to be included.
#' @param colNames A logical specifying if the column names are to be included.
#' @param upperLeft A character string specifying the entry at the upper left corner.
#' @param class A character string specifying the HTML global class attribute.
#'
#' @examples
#' myDat <- GE.read.csv(system.file("extdata", "F2maize_pheno.csv", package = "RAP"),
#'                      env = "env!", genotype = "genotype!", trait = "yld")
#' myTD <- createTD(data = myDat, genotype = "genotype!", env = "env!")
#' stats <- summary(object = myTD, trait = "yld", all = TRUE)
#' RAP.htmlTable(stats)
#'
#' @export

RAP.htmlTable <- function(x,
                          rowNames = TRUE,
                          colNames = TRUE,
                          upperLeft = "",
                          class = "GenTable") {
  x <- as.matrix(x)
  d0 <- dim(x)
  d1 <- d0[1]
  d2 <- d0[2]
  rNames <- rownames(x)
  cNames <- colnames(x)
  # Create a string to store html table
  hTable <- paste0("<TABLE CLASS=", class, ">\n")
  # Write header of the table
  if (colNames && !is.null(cNames)) {
    if (length(cNames) != d2) {
      stop("The number of column names does not equal the number of columns.")
    }
    if (rowNames && !is.null(rNames)) {
      if (length(rNames) != d1) {
        stop("The number of row names does not equal the number of rows.")
      }
      for (i in 0:d2) {
        if (i == 0) {
          hTable <- paste0(hTable, "<TR><TD>", upperLeft, "</TD>\n")
        } else if (i == d2) {
          hTable <- paste0(hTable, "<TD>", cNames[i], "</TD></TR>\n")
        } else {
          hTable <- paste0(hTable, "<TD>", cNames[i], "</TD>\n")
        }
      }
    } else {
      for (i in 1:d2) {
        if (i == 1) {
          hTable <- paste0(hTable, "<TR><TD>", cNames[i], "</TD>\n")
        } else if (i == d2) {
          hTable <- paste0(hTable, "<TD>", cNames[i], "</TD></TR>\n")
        } else {
          hTable <- paste0(hTable, "<TD>", cNames[i], "</TD>\n")
        }
      }
    }
  }
  # Write contents of the table
  if (rowNames && !is.null(rNames)) {
    for (i in 1:d1) {
      for (j in 1:d2) {
        if (j == 1) {
          hTable <- paste0(hTable, "<TR><TD>", rNames[i], "</TD><TD>&nbsp;", x[i, j], "</TD>\n")
        }else if (j==d2) {
          hTable <- paste0(hTable, "<TD>&nbsp;", x[i, j], "</TD></TR>\n")
        } else {
          hTable <- paste0(hTable, "<TD>&nbsp;", x[i, j], "</TD>\n")
        }
      }
    }
  } else {
    for (i in 1:d1) {
      for (j in 1:d2) {
        if (j == 1){
          hTable <- paste0(hTable,"<TR><TD>&nbsp;", x[i, j],"</TD>\n")
        } else if (j == d2) {
          hTable <- paste0(hTable,"<TD>&nbsp;",x[i,j],"</TD></TR>\n")
        } else{
          hTable <- paste0(hTable,"<TD>&nbsp;",x[i,j],"</TD>\n")
        }
      }
    }
  }
  return(paste(hTable, "</TABLE>\n"))
}
