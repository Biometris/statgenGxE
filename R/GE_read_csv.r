#' Data Input
#'
#' This function reads a comma-separated values file for a GxE analysis.
#'
#' @param csvfile The name of the comma-separated values file which the data are to be read from.
#' @param env A character string specifying the column name of environments.
#' @param genotype A character string specifying the column name of the genotypes.
#' @param trait A character string specifying a trait name.
#' @param megaEnv A character vector specifying mega environment factor names. Default, \code{NULL}.
#' @param sep The field separator character. Values on each line of the file are separated by this character.
#'   If sep = "" (the default for read.table) the separator is `white space', that is one or more spaces, tabs, newlines or carriage returns.
#' @param quote The set of quoting characters. To disable quoting altogether, use quote = "". See \code{\link{scan}} for the behaviour on quotes embedded in quotes.
#'   Quoting is only considered for columns read as character, which is all of them unless colClasses is specified.
#' @param dec The character used in the file for decimal points.
#' @param fill A logical. If TRUE, when the rows have unequal length blank fields are implicitly added. See details in \code{\link{read.table}}.
#' @param comment.char A character vector of length one containing a single character or an empty string. Use "" to turn off the interpretation of comments altogether.
#' @param envSelect A character vector of unique names to be selected from the environment ID(s). If \code{NA}, all the rows are included.
#' @param na.strings A character vector of strings which are to be interpreted as \code{NA} values.
#'   Blank fields are also considered to be missing values in logical, integer, numeric and complex fields.
#' @param check.names A logical. If \code{TRUE} then the names of the variables in the data frame are checked to ensure that they are syntactically valid variable names.
#'   If necessary they are adjusted (by \code{\link{make.names}}) so that they are, and also to ensure that there are no duplicates.
#' @param ... Further arguments to be passed to read.table.
#' @return A data frame (\code{\link{data.frame}}) containing a representation of (a sebset of)the data in the file.
#'
#' @seealso \code{\link{read.csv}} which this function wraps.
#' @export
#' @note a header is needed in the csv file which the data are to be read from.
#' @examples
#' mydat <- GE.read.csv(file.path(path.package("RAP"),"F2maize_pheno.csv"),
#'                      env="env!", genotype="genotype!", trait="yld")
#' str(mydat)
#' 
#' @export
GE.read.csv <- function(csvfile, env, genotype, trait, megaEnv=NULL, sep = ",",
            quote="\"", dec=".", fill = TRUE, comment.char="", envSelect=NA, na.strings = "NA", check.names=FALSE, ...){

    # read the header lines as char strings
    hl=readLines(csvfile, 1)

    # split headers up by 'sep'
    hl=unlist(strsplit(hl,sep))
    hl=gsub(quote,"",hl)
    nhl=length(hl)
    if(nhl!=length(unique(hl)))
      stop("There are some duplicates of the names on the header.\n")
    cls=rep(NA,nhl)

    # the variables to be factorized
    cls[env==hl]="factor"
    cls[genotype==hl]="factor"
    # the variables to be numeric
    cls[hl%in%trait]="numeric"
    if(!is.null(megaEnv)){
      cls[hl%in%megaEnv] <- "factor"
    }

    # use read.csv to read the data set
    mydat<-read.csv(csvfile, header = TRUE, sep = sep, quote=quote, dec=dec,
         fill = fill, comment.char=comment.char,colClasses=cls, na.strings = na.strings, check.names=check.names, ...)

    if (!(env %in% names(mydat))) stop(env, " not found!\n")
    if (!(genotype %in% names(mydat))) stop(genotype, " not found!\n")
    for (trt in trait)
      if (!(trt %in% names(mydat))) stop(trt, " not found!\n")

    if (!all(is.na(envSelect))){
        mydat=mydat[mydat[[env]]%in%envSelect,]
        mydat=droplevels(mydat)
        # check if the number of rows is correct
        if(length(unique(mydat[[env]]))!= length(envSelect))
          warning("CHECK: Either some selected row(s) are not read in\n or the selected row name(s) do not match the enviroment ID(s)!")
    }

    # validating data for GxE analysis
    ngeno <- nlevels(mydat[[genotype]])
    nenv  <- nlevels(mydat[[env]])
    ntrait <- nrow(mydat)

    # check if the supplied data contains the genotype by environment means
    if (ntrait != ngeno * nenv)
      message("There is not a trait mean per genotype per enviroment")
    # check if any missing value in envrironment and genotype
    if (any(is.na(mydat[,genotype])))
      message("some missing value(s) found in ", genotype, " column")
    if (any(is.na(mydat[,env])))
      message("some missing value(s) found in ",      env, " column")

    # makes sure that a regular 2-d grid is produced.
    if ((ntrait != ngeno * nenv)|any(is.na(mydat[,genotype]))|any(is.na(mydat[,env]))){
      temptab <- tapply(1:ntrait, mydat[,c(genotype, env)], identity)
      ind <- which(is.na(temptab), arr.ind=T)
      if (length(ind)>0){
        addgenonames <- rownames(temptab)[ind[,1]]
        addenvnames  <- colnames(temptab)[ind[,2]]
      }
      if (any(is.na(mydat[,genotype]))|any(is.na(mydat[,env]))){
        delpos <- c(which(is.na(mydat[,genotype])),which(is.na(mydat[,env])))
        delpos <- unique(delpos)
        mydat <- mydat[-delpos,]
        rownames(mydat) <- 1:nrow(mydat)
      }
      if (length(ind)>0){
        nr <- nrow(mydat)
        for (i in 1:length(addgenonames)){
          mydat[(nr+i),genotype] <- addgenonames[i]
          mydat[(nr+i),env] <- addenvnames[i]
        }
      }
      mydat=droplevels(mydat)
      message("expands/reduces vectors onto a regular 2-d grid")
    }

    # removing any environment or genotype with completely missing values
    if (any(is.na(mydat[,trait]))){
      test1 <- sapply(1:length(trait), function(y) tapply(mydat[[trait[y]]], mydat[,genotype], function(x) all(is.na(x))))
      test1 <- apply(test1, 1, all)
      if (any(test1)){
        delnames <- names(test1)[test1]
        delpositions <- mydat[[genotype]] %in% delnames
        mydat=mydat[!delpositions,]
        mydat=droplevels(mydat)
        message("Genotype(s) ", paste(delnames, collapse=",") ," with all completely missing values are removed.")
      }
      test2 <- sapply(1:length(trait), function(y) tapply(mydat[[trait[y]]], mydat[,env], function(x) all(is.na(x))))
      test2 <- apply(test2, 1, all)
      if (any(test2)){
        delnames2 <- names(test2)[test2]
        delpositions2 <- mydat[[env]] %in% delnames
        mydat=mydat[!delpositions2,]
        mydat=droplevels(mydat)
        message("Environment(s) ",paste(delnames2, collapse=",") ," with all completely missing values are removed.")
      }
    }
    # return data in the specific order
    mydat <- mydat[order(mydat[[genotype]],mydat[[env]]),]
    rownames(mydat) <- 1:nrow(mydat)
    mydat
}
