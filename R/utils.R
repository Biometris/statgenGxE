#' Custom tryCatch to return result, errors and warnings.
#' Copied from http://stackoverflow.com/a/24569739/2271856.
#'
#' @keywords internal
tryCatchExt <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }), warning = function(w) {
      warn <<- c(warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(value = value, warning = warn, error = err)
}

#' Helper function for suppressing a single warning message.
#' @keywords internal
supprWarn <- function(expression,
                      message) {
  withCallingHandlers(expression,
                      warning = function(w) {
                        if (grepl(message, w$message)) {
                          invokeRestart("muffleWarning")
                        }
                      })
}

#' Helper function for checking whether error message about 1% change on
#' last iteration for asreml is worth mentioning as a warning.
#' If the corresponding parameter is close to zero and then changes off 1%
#' or more can be expected an are ok.
#' @keywords internal
chkLastIter <- function(model) {
  wrnMsg <- paste("At least one parameter changed by more than 1%",
                  "on the last iteration")
  if (any(grepl(pattern = wrnMsg, x = model$warning))) {
    ## EXtract monitor df from model object.
    mon <- model$value$monitor
    ## Extract values for parameters for last 2 iterations.
    ## First 3 rows give general model info. Last col a summary.
    lastIt <- mon[-(1:3), c(ncol(mon) - 2, ncol(mon) - 1)]
    ## Compute change of parameters in last iteration.
    change <- abs((lastIt[, 2] - lastIt[, 1]) / lastIt[, 1]) * 100
    ## Suppress waning if the change was less than 5% or the param value less
    ## than 0.1.
    if (all(change <= 5) || all(lastIt[change > 5, 1] < 0.1)) {
      model$warning <- model$warning[!grepl(pattern = wrnMsg,
                                            x = model$warning)]
    }
  }
  return(model)
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
  wrnMsg <- "Insufficient workspace - (reset workspace or pworkspace arguments)"
  ## Create tempfile to suppress asreml output messages.
  tmp <- tempfile()
  sink(tmp)
  ## Predict using default settings, i.e. pworkspace = 8e6
  modelP <- tryCatchExt(predict(model, classify = classify,
                                vcov = vcov, associate = associate,
                                data = TD, maxiter = 20, ...))
  pWorkSpace <- 8e6
  ## While there is a warning, increase pWorkSpace and predict again.
  while (!is.null(modelP$warning) &&
         any(grepl(pattern = wrnMsg, x = modelP$warning))
         && pWorkSpace < 160e6) {
    pWorkSpace <- pWorkSpace + 8e6
    modelP <- tryCatchExt(predict(model, classify = classify,
                                  vcov = vcov, associate = associate, data = TD,
                                  maxiter = 20, pworkspace = pWorkSpace, ...))
  }
  sink()
  unlink(tmp)
  if (!is.null(modelP$warning) && !all(grepl(pattern = wrnMsg,
                                             x = modelP$warning))) {
    modelP <- chkLastIter(modelP)
    if (length(modelP$warning) != 0) {
      warning(modelP$warning, call. = FALSE)
    }
  }
  if ((length(modelP$warning) == 0 ||
       !all(grepl(pattern = wrnMsg, x = modelP$warning))) &&
      is.null(modelP$error)) {
    return(modelP$value)
  } else {
    stop(paste("Error in asreml when running predict. Asreml message:\n",
               modelP$error, "\n",
               modelP$warning, "\n"), call. = FALSE)
  }
}

#' Helper function for computing the standard error of the variance.
#'
#' @keywords internal
seVar <- function(x,
                  na.rm = FALSE) {
  if (inherits(x, c("matrix", "data.frame"))) {
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
    se <- sqrt((n * (m4 - 4 * m1 * m3 + 6 * m1 ^ 2 * m2 - 3 * m1 ^ 4) /
                  (n - 1) - (n * (m2 - m1 ^ 2) / (n - 1)) ^ 2) / n)
  } else {
    se <- seVar(x = as.vector(x), na.rm = na.rm)
  }
  return(se)
}

#' Helper function for computing the skewness.
#' This and following formulas taken from
#' https://brownmath.com/stat/shape.htm#Normal.
#'
#' @keywords internal
skewness <- function(x,
                     na.rm = FALSE) {
  if (inherits(x, c("matrix", "data.frame"))) {
    skw <- apply(X = x, MARGIN = 2, FUN = skewness, na.rm = na.rm)
  } else if (is.vector(x)) {
    if (na.rm) {
      x <- x[!is.na(x)]
    }
    n <- length(x)
    skw <- (sum((x - mean(x)) ^ 3) / n) / (sum((x - mean(x)) ^ 2) / n) ^ (3 / 2)
  } else {
    skw <- skewness(x = as.vector(x), na.rm = na.rm)
  }
  return(skw)
}

#' Helper function for computing the standard error of the skewness.
#'
#' @keywords internal
seSkewness <- function(n) {
  if (n <= 2) {
    warning(paste("For n less than 2 the standard error of skewness cannot be",
                  "calculated"), call. = FALSE)
    return(NA)
  }
  return(sqrt((6 * n * (n - 1)) / ((n - 2) * (n + 1) * (n + 3))))
}

#' Helper function for computing kurtosis.
#' Rescaled by subtracting 3 from the result to give the normal distribution
#' a kurtosis of 0, so basically the excess kurtosis.
#'
#' @keywords internal
kurtosis <- function(x,
                     na.rm = FALSE) {
  if (inherits(x, c("matrix", "data.frame"))) {
    kurt <- apply(X = x, MARGIN = 2, FUN = kurtosis, na.rm = na.rm)
  } else if (is.vector(x)) {
    if (na.rm) {
      x <- x[!is.na(x)]
    }
    n <- length(x)
    kurt <- n * sum((x - mean(x)) ^ 4) / (sum((x - mean(x)) ^ 2) ^ 2) - 3
  } else {
    kurt <- kurtosis(x = as.vector(x), na.rm = na.rm)
  }
  return(kurt)
}

#' Helper function for computing the standard error of the kurtosis.
#'
#' @keywords internal
seKurtosis <- function(n) {
  if (n <= 3) {
    warning(paste("For n less than 2 the standard error of kurtosis cannot be",
                  "calculated"), call. = FALSE)
    return(NA)
  }
  return(sqrt((24 * n * (n - 1) ^ 2) / ((n - 2) * (n - 3) * (n + 3) * (n + 5))))
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
           ## Construct shell command for changing directory.
           ## cd /d is used instead of cd to account for changing drives on windows.
           ## Note that here dirname(outfile) is needed instead of outDir.
           cmdDir <- paste0("cd /d ", dirname(outfile))
           shell(cmd = paste(cmdDir, "&", cmdRun1, "> nul 2>&1"))
           shell(cmd = paste(cmdDir, "&", cmdRun2, "> nul"))
         }, linux = {
           ## Construct shell command for changing directory.
           cmdDir <- paste("cd", outDir)
           system(command = paste(cmdDir, ";", cmdRun1, "> /dev/null 2>&1"))
           system(command = paste(cmdDir, ";", cmdRun2, "> /dev/null"))
         }, darwin = {
           ## Construct shell command for changing directory.
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
    regmatches(x = unlist(strsplit(cp, "@"))[2],
               gregexpr("[[:digit:]]+\\.*[[:digit:]]",
                        unlist(strsplit(cp, "@"))[2]))[[1]]
  })
  posExt <- sapply(X = chrPos, function(cp) {
    sub(pattern = "\\d*(\\.\\d)", replacement = "",
        x = unlist(strsplit(cp, "@"))[2])
  })
  chrPosDf <- data.frame(chr, pos, posExt, stringsAsFactors = FALSE)
  mapTot <- qtl::pull.map(cross, as.table = TRUE)[, 1:2]
  mapTot$mrkNames <- rownames(mapTot)
  ## If there is an X chromosome, position is named pos.female. Rename for
  ## easier merging. Pos.male is ignored.
  colnames(mapTot)[2] <- "pos"
  mapTot$pos <- round(mapTot$pos, digits = 1)
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



