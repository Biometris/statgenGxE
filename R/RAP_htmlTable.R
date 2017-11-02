#' Generates an HTML table
#'
#' Write a matrix or data frame as an HTML table.
#'
#' @param x A matrix or data frame object.
#' @param row.names A logical specifying if the row names are to be included.
#' @param col.names A logical specifying if the column names are to be included.
#' @param upper.left A character string specifying the entry at the upper left corner.
#' @param class A character string specifying the HTML global class attribute.
#' @examples
#' mydat <- ST.read.csv(file.path(path.package("RAP"),"SB_yield.csv"),
#'                      factor.names=c("Env","Genotype","Rep","Subblock","Row","Column"),
#'                      trait.names="yield", env ="Env", rowSelect="HEAT06",
#'                      colSelect=c("Env","Genotype","Rep","Row","Column","yield"))
#' stats <- ST.summary.trait(data=mydat, trait="yield", all=TRUE)
#' RAP.htmlTable(stats)
#'
#' @export


RAP.htmlTable <- function(x, row.names=TRUE, col.names=TRUE, upper.left="", class="GenTable"){

  x <- as.matrix(x)
  d0 <- dim(x)
  d1 <- d0[1]
  d2 <- d0[2]
  rnames <- rownames(x)
  cnames <- colnames(x)

  # Create a string to store html table
  htable <- paste0("<TABLE CLASS=", class, ">\n")
  # Write header of the table
  if (isTRUE(col.names) && !is.null(cnames)){
    if (length(cnames)!=d2)
      stop("The number of column names does not equal the number of columns.")
    if (isTRUE(row.names) && !is.null(rnames)){
      if (length(rnames)!=d1)
        stop("The number of row names does not equal the number of rows.")
      for (i in 0:d2){
        if (i==0){
          htable <- paste(htable,"<TR><TD>", upper.left,"</TD>\n",sep="")
        }else{
          if (i==d2){
            htable <- paste(htable,"<TD>",cnames[i],"</TD></TR>\n",sep="")
          }else{
            htable <- paste(htable,"<TD>",cnames[i],"</TD>\n",sep="")
          }
        }
      }
    }else{
      for (i in 1:d2){
        if (i==1){
          htable <- paste(htable,"<TR><TD>",cnames[i],"</TD>\n",sep="")
        }else{
          if (i==d2){
            htable <- paste(htable,"<TD>",cnames[i],"</TD></TR>\n",sep="")
          }else{
            htable <- paste(htable,"<TD>",cnames[i],"</TD>\n",sep="")
          }
        }
      }
    }
  }

  # Write contents of the table
  if (isTRUE(row.names) && !is.null(rnames)){
    for (i in 1:d1){
      for (j in 1:d2){
        if (j==1){
          htable <- paste(htable,"<TR><TD>", rnames[i],"</TD><TD>&nbsp;",x[i,j],"</TD>\n",sep="")
        }else{
          if (j==d2){
            htable <- paste(htable,"<TD>&nbsp;",x[i,j],"</TD></TR>\n",sep="")
          }else{
            htable <- paste(htable,"<TD>&nbsp;",x[i,j],"</TD>\n",sep="")
          }
        }
      }
    }
  }else{
    for (i in 1:d1){
      for (j in 1:d2){
        if (j==1){
          htable <- paste(htable,"<TR><TD>&nbsp;",x[i,j],"</TD>\n",sep="")
        }else{
          if (j==d2){
            htable <- paste(htable,"<TD>&nbsp;",x[i,j],"</TD></TR>\n",sep="")
          }else{
            htable <- paste(htable,"<TD>&nbsp;",x[i,j],"</TD>\n",sep="")
          }
        }
      }
    }
  }
  paste(htable, "</TABLE>\n",sep="")
}
