#' Custom tryCatch to return result, errors and warnings.
#' Copied from http://stackoverflow.com/a/24569739/2271856.
#'
#' @keywords internal
tryCatchExt <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error = function(e) {
      err <<- e
      NULL
    }), warning = function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(value = value, warning = warn, error = err)
}

#' Extended version of asreml.predict
#'
#' Asreml has a bug that may throw a warning message:
#' Abnormal termination
#' Insufficient workspace - (reset workspace or pworkspace arguments)
#' This may be avoided by increasing pworkspace, but this doesn't
#' always work.
#' If this happens pworkspace is increased in 'small' steps.
#'
#' @keywords internal
predictAsreml <- function(model,
                          classify = "genotype",
                          associate = as.formula("~ NULL"),
                          vcov = TRUE,
                          TD,
                          ...) {
  ## Create tempfile to suppress asreml output messages.
  tmp <- tempfile()
  sink(tmp)
  ## Predict using default settings, i.e. pworkspace = 8e6
  modelP <- tryCatchExt(predict(model, classify = classify,
                                vcov = vcov, associate = associate,
                                data = TD, ...))
  pWorkSpace <- 8e6
  ## While there is a warning, increase pWorkSpace and predict again.
  while (!is.null(modelP$warning) && pWorkSpace < 160e6) {
    pWorkSpace <- pWorkSpace + 8e6
    modelP <- tryCatchExt(predict(model, classify = classify,
                                  vcov = vcov, associate = associate, data = TD,
                                  pworkspace = pWorkSpace, ...))
  }
  sink()
  unlink(tmp)
  if (is.null(modelP$warning) && is.null(modelP$error)) {
    return(modelP$value)
  } else {
    stop(paste("Error in asreml when running predict. Asreml message:\n",
               modelP$error$message, "\n",
               modelP$warning$message, "\n"), call. = FALSE)
  }
}

#' @keywords internal
checkCols <- function(cols,
                      data) {
  cCheck <- !cols %in% colnames(data)
  if (any(cCheck)) {
    stop(paste0("Error in ", deparse(sys.call(-1)), ":\n\t",
                "The following colomns are not in ",
                deparse(substitute(data, env = parent.frame(2))),
                ": ", paste(cols[cCheck], collapse = ", "), "\n"),
         call. = FALSE)
  }
}

#' @keywords internal
isValidVariableName <- function(x,
                                allowReserved = TRUE,
                                unique = FALSE) {
  # Author: Richie Cotton
  # http://4dpiecharts.com/tag/regex/
  ok <- rep.int(x = TRUE, times = length(x))
  #is name too long?
  #max_name_length <- if(getRversion() < "2.13.0") 256L else 10000L
  #ok[nchar(x) > max_name_length] <- FALSE
  #is it a reserved variable, i.e.
  #an ellipsis or two dots then a number?
  if (!allowReserved) {
    ok[x == "..."] <- FALSE
    ok[grepl(pattern = "^\\.{2}[[:digit:]]+$", x = x)] <- FALSE
  }
  #are names valid (and maybe unique)
  ok[x != make.names(x, unique = unique)] <- FALSE
  return(ok)
}

#' @keywords internal
char2numeric <- function(x,
                         dec = ".") {
  # This function is to convert a vector of characters to numeric with specified dec
  if (dec == ".") {
    y <- as.numeric(x)
  } else {
    y <- sapply(X = strsplit(x, dec, fixed = TRUE), FUN = function(z) {
      as.numeric(paste(z, collapse = "."))
    })
  }
  return(y)
}

#' @keywords internal
seVar <- function(x, na.rm = FALSE)
{
  if (is.matrix(x)) {
    se <- apply(X = x, MARGIN = 2, FUN = seVar, na.rm = na.rm)
  } else if (is.vector(x)) {
    if (na.rm) {
      x <- x[!is.na(x)]
    }
    n <- length(x)
    m1 <- sum(x) / n
    m2 <- sum(x ^ 2) / n
    m3 <- sum(x ^ 3) / n
    m4 <- sum(x ^ 4) / n
    se <- sqrt((n * (m4 - 4 * m1 * m3 + 6 * m1 ^ 2 * m2 -
                       3 * m1 ^ 4) / (n - 1) - (n * (m2 - m1 ^ 2) /
                                                  (n - 1)) ^ 2) / n)
  } else if (is.data.frame(x)) {
    se <- sapply(X = x, FUN = seVar, na.rm = na.rm)
  } else {
    se <- seVar(x = as.vector(x), na.rm = na.rm)
  }
  return(se)
}

#' @keywords internal
skewness <- function(x,
                     na.rm = FALSE) {
  if (is.matrix(x)) {
    skw <- apply(X = x, MARGIN = 2, FUN = skewness, na.rm = na.rm)
  } else if (is.vector(x)) {
    if (na.rm) {
      x <- x[!is.na(x)]
    }
    n <- length(x)
    m1 <- sum(x) / n
    m2 <- sum(x ^ 2) / n
    m3 <- sum(x ^ 3) / n
    skw <- (m3 - 3 * m1 * m2 + 2 * m1 ^ 3) / (m2 - m1 ^ 2) ^ (3 / 2)
  } else if (is.data.frame(x)) {
    skw <- sapply(X = x, FUN = skewness, na.rm = na.rm)
  } else {
    skw <- skewness(x = as.vector(x), na.rm = na.rm)
  }
  return(skw)
}

#' @keywords internal
seSkewness <- function(n) {
  return(sqrt((6 * n * (n - 1)) / ((n - 1) * (n + 1) * (n + 3))))
}

#' @keywords internal
kurtosis <- function(x,
                     na.rm = FALSE) {
  if (is.matrix(x)) {
    kurt <- apply(X = x, MARGIN = 2, FUN = kurtosis, na.rm = na.rm)
  } else if (is.vector(x)) {
    if (na.rm) {
      x <- x[!is.na(x)]
    }
    n <- length(x)
    m1 <- sum(x) / n
    m2 <- sum(x ^ 2) / n
    m3 <- sum(x ^ 3) / n
    m4 <- sum(x ^ 4) / n
    kurt <- (m4 - 4 * m1 * m3 + 6 * m1 ^ 2 * m2 - 3 * m1 ^ 4) /
      (m2 - m1 * m1) ^ 2 - 3
  } else if (is.data.frame(x)) {
    kurt <- sapply(X = x, FUN = kurtosis, na.rm = na.rm)
  } else {
    kurt <- kurtosis(x = as.vector(x), na.rm = na.rm)
  }
  return(kurt)
}

#' @keywords internal
seKurtosis <- function(n) {
  return(sqrt((24 * n * (n - 1) ^ 2) / ((n - 2) * (n - 3) * (n + 5) * (n + 3))))
}

#' Base method for creating a report
#'
#' Base method for creating a .pdf and .tex report from an \code{R} object
#'
#' @param x An \code{R} object
#' @param ... Further arguments to be passed on to specific report functions.
#'
#' @export
report <- function(x,
                   ...) {
  UseMethod("report")
}

#' Helper function for creating the actual report
#'
#' @keywords internal
createReport <- function(x,
                         reportName,
                         outfile,
                         ...) {
  ## Check provided outfile
  if (!is.null(outfile)) {
    if (!is.character(outfile) || length(outfile) > 1 ||
        tools::file_ext(outfile) != "pdf") {
      stop("invalid output filename provided.\n")
    }
    ## Since latex cannot handle spaces in figure paths knitr converts those
    ## to pathnames with _. To prevent this spaces are not allowed.
    if (grepl(pattern = " ", x = outfile)) {
      stop("outfile path cannot contain spaces. Provide a path without spaces or
           a relative path")
    }
  } else {
    ## Create a generic output filenane from the name of the report and
    ## the current date/time. The file will be placed in the current
    ## working directory.
    timeStamp <- format(Sys.time(), "%Y%m%d%H%M%S")
    outfile <- paste0("./" , substring(reportName, first = 1,
                                       last = nchar(reportName) - 4),
                      "_", timeStamp, ".pdf")
  }
  ## Extract output directory from outfile.
  outDir <- dirname(outfile)
  ## If output directory doesn't exist, create it.
  if (!dir.exists(outDir)) {
    dir.create(outDir, recursive = TRUE)
  }
  ## When the output need to be written to a top level directory on windows
  ## there may be an extra / at the end of the filename.
  ## This is removed here so file.path works properly further on.
  if (tolower(Sys.info()[["sysname"]]) == "windows") {
    if (substring(text = outDir,
                  first = nchar(outDir)) == .Platform$file.sep) {
      outDir <- substring(text = outDir, first = 1,
                          last = nchar(outDir) - 1)
    }
  }
  ## Extract the name of the outputfile, so without path and extension.
  outBase <- substring(basename(outfile), first = 1,
                       last = nchar(basename(outfile)) - 3)
  ## Construct the output name of the .tex file
  outTex <- file.path(outDir, paste0(outBase, "tex"))
  ## Get the report file from the directory where the package is installed.
  reportFile <- system.file("reports", reportName, package = "RAP")
  ## Save knitr options for reset when exiting function.
  knitrOptsOrig <- knitr::opts_chunk$get()
  on.exit(knitr::opts_chunk$set(knitrOptsOrig))
  ## Run knitr with chunk options set to produce proper ppt.
  figPrefix <- paste0(format(Sys.time(), "%m%d%y%H%M%S"), "-")
  knitr::opts_chunk$set(fig.show = "hold",
                        fig.path = file.path(outDir, "figures", figPrefix),
                        fig.process = function(x) {
                          paste0("./figures/", basename(x))
                        })
  knitr::knit(input = reportFile, output = outTex, quiet = TRUE)
  ## Construct shell commands for calling pdf latex.
  ## First only draftmode for speed.
  cmdRun1 <- paste0("pdflatex -interaction=nonstopmode -draftmode ",
                    outBase, "tex")
  cmdRun2 <- paste0("pdflatex -interaction=nonstopmode ",
                    outBase, "tex")
  ## Run shell commands. System doens't work for windows.
  ## Two runs needed to get references right.
  switch(tolower(Sys.info()[["sysname"]]),
         windows = {
           ## Construct shell command for changing directory
           ## cd /d is used instead of cd to account for changing drives on windows.
           ## Note that here dirname(outfile) is needed instead of outDir.
           cmdDir <- paste0("cd /d ", dirname(outfile))
           shell(cmd = paste(cmdDir, "&", cmdRun1, "> nul 2>&1"))
           shell(cmd = paste(cmdDir, "&", cmdRun2, "> nul"))
         }, linux = {
           ## Construct shell command for changing directory
           cmdDir <- paste("cd", outDir)
           system(command = paste(cmdDir, ";", cmdRun1, "> /dev/null 2>&1"))
           system(command = paste(cmdDir, ";", cmdRun2, "> /dev/null"))
         }, darwin = {
           ## Construct shell command for changing directory
           cmdDir <- paste("cd", outDir)
           system(command = paste(cmdDir, ";", cmdRun1, "> /dev/null 2>&1"))
           system(command = paste(cmdDir, ";", cmdRun2, "> /dev/null"))
         })
  ## Remove extra files generated by pdflatex.
  for (extension in c("aux", "log", "out", "toc", "xwm")) {
    unlink(file.path(outDir, paste0(outBase, extension)))
  }
  ## A .csv file might be created in the package directory.
  ## This file is moved to the proper output location.
  if (file.exists(system.file("reports", paste0(outBase, "csv"),
                              package = "RAP"))) {
    file.copy(system.file("reports", paste0(outBase, "csv"), package = "RAP"),
              file.path(outDir, paste0(outBase, "csv")), overwrite = TRUE)
    unlink(system.file("reports", paste0(outBase, "csv"), package = "RAP"))
  }
}

#' @keywords internal
qtlPosToName <- function(chrPos, cross) {
  chr <- sapply(X = chrPos, function(cp) {
    unlist(strsplit(cp, "@"))[1]
  })
  pos <- sapply(X = chrPos, function(cp) {
    as.numeric(unlist(strsplit(unlist(strsplit(cp, "@"))[2], "[[:alpha:]]"))[1])
  })
  posExt <- sapply(X = chrPos, function(cp) {
    unlist(strsplit(unlist(strsplit(cp, "@"))[2],
                    "(?=[A-Za-z])(?<=[0-9])", perl = TRUE))[2]
  })
  chrPosDf <- data.frame(chr, pos, posExt, stringsAsFactors = FALSE)
  mapTot <- qtl::pull.map(cross, as.table = TRUE)
  mapTot$mrkNames <- rownames(mapTot)
  chrPosMap <- merge(chrPosDf, mapTot, by = c("chr", "pos"), all.x = TRUE)
  chrPosMap[is.na(chrPosMap$mrkNames), "mrkNames"] <-
    paste0("c", chrPosMap[is.na(chrPosMap$mrkNames), "chr"], ".loc",
           chrPosMap[is.na(chrPosMap$mrkNames), "pos"])
  return(list(chrNames = chrPosMap$mrkNames, ext = chrPosMap$posExt))
}

#' Function for plotting a correlation (or covariance) matrix.
#
#' @keywords internal
plotCorMat <- function(corMat, main = "") {
  meltedCorMat <- reshape2::melt(corMat)
  ggplot2::ggplot(data = meltedCorMat, ggplot2::aes_string("Var1", "Var2",
                                                           fill = "value")) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                                  na.value = "grey", midpoint = 0,
                                  limit = c(-1, 1), space = "Lab") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1,
                                                       size = 10, hjust = 1)) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::ggtitle(main) + ggplot2::xlab("") + ggplot2::ylab("") +
    ggplot2::coord_fixed()
}



