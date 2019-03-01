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
    change <- ifelse(lastIt[, 1] == 0, 0, abs((lastIt[, 2] - lastIt[, 1]) /
                                                lastIt[, 1]) * 100)
    ## Suppress warning if the change was less than 5% or the param value less
    ## than 0.1.
    if (all(change <= 5) || all(lastIt[change > 5, 1] < 0.1)) {
      model$warning <- model$warning[!grepl(pattern = wrnMsg,
                                            x = model$warning)]
    }
  }
  return(model)
}

#' Helper function for converting certain asreml warnings to errors.
#' @keywords internal
wrnToErr <- function(model) {
  wrns <- c("Abnormal termination", "returning -Inf")
  for (wrn in wrns) {
    if (any(grepl(pattern = wrn, x = model$warning))) {
      ## Remove from warnings and add to errors
      model$error <- c(model$error, model$warning[grepl(pattern = wrn,
                                                        x = model$warning)])
      model$warning <- model$warning[!grepl(pattern = wrn,
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
  wrnMsg <- "reset workspace or pworkspace arguments"
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
#' Base method for creating a .pdf and .tex report from an \code{R} object.
#'
#' @param x An \code{R} object
#' @param ... Further arguments to be passed on to specific report functions.
#'
#' @seealso \code{\link{report.SSA}}, \code{\link{report.varComp}},
#' \code{\link{report.AMMI}}, \code{\link{report.FW}},
#' \code{\link{report.stability}}, \code{\link{report.cross}},
#' \code{\link{report.QTLDet}}, \code{\link{report.multiQTL}}
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
  cmdRun1 <- paste0(Sys.which("pdflatex"), " -interaction=nonstopmode -draftmode ",
                    outBase, "tex")
  cmdRun2 <- paste0(Sys.which("pdflatex"), " -interaction=nonstopmode ",
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

#' Function for extracting the name of a qtl combining the peaks from a QTLDet
#' object and the position of the qtl on the chromosome.
#'
#' @keywords internal
qtlPosToName <- function(chrPos,
                         peaks) {
  ## chromosome positions are given in chr@position + .a, .d or nothing.
  ## Extract chromosome.
  chr <- sapply(X = chrPos, function(cp) {
    unlist(strsplit(cp, "@"))[1]
  })
  ## Extract chromosome position in two steps. First take the part behind @.
  ## Then remove everything behind the dot.
  pos <- as.numeric(sapply(X = chrPos, function(cp) {
    regmatches(x = unlist(strsplit(cp, "@"))[2],
               gregexpr("[[:digit:]]+\\.*[[:digit:]]",
                        unlist(strsplit(cp, "@"))[2]))[[1]]
  }))
  ## Extract the extesion, i.e. the bit after the last dot.
  posExt <- sapply(X = chrPos, function(cp) {
    sub(pattern = "\\d*(\\.\\d)", replacement = "",
        x = unlist(strsplit(cp, "@"))[2])
  })
  chrPosDf <- data.frame(chr, pos, posExt, stringsAsFactors = FALSE)
  #mapTot <- qtl::pull.map(cross, as.table = TRUE)[, 1:2]
  mapTot <- peaks[ , c("chr", "pos")]
  mapTot$mrkNames <- rownames(mapTot)
  ## QTLDetection rounds positions to 1 digit in the names.
  ## For a proper match position in the peaks has to be rounded to 1
  ## digit as well.
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
plotCorMat <- function(varMat,
                       main = "") {
  ## Melt variance and correlation matrices to get proper shape for ggplot.
  meltedCorMat <- reshape2::melt(cov2cor(varMat))
  meltedvarMat <- reshape2::melt(varMat)
  ## Select bottom triangle for correlations and top for variances.
  meltedCorMatLow <- meltedCorMat[as.numeric(meltedCorMat$Var1) >
                                    as.numeric(meltedCorMat$Var2), ]
  meltedvarMatUp <- meltedvarMat[as.numeric(meltedvarMat$Var1) <=
                                   as.numeric(meltedvarMat$Var2), ]
  ## Round values for nicer display
  meltedvarMatUp$value <- round(meltedvarMatUp$value)
  ggplot2::ggplot(data = meltedCorMatLow,
                  ggplot2::aes_string("Var1", "Var2", fill = "value")) +
    ggplot2::geom_tile(color = "white") +
    ## Create a gradient scale.
    ggplot2::scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                                  na.value = "grey", limit = c(-1, 1)) +
    ggplot2::geom_text(data = meltedvarMatUp,
                       ggplot2::aes_string(label = "value", size = "value")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1,
                                                       size = 10, hjust = 1)) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ## Remove grid behind text output.
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::ggtitle(main) + ggplot2::xlab("") + ggplot2::ylab("") +
    ggplot2::labs(fill = "correlation", size = "variance") +
    ## Fix coordinates to get a square sized plot.
    ggplot2::coord_fixed()
}

#' Function for escaping special LaTeX characters
#'
#' Taken from knitr package. Copied since it is an internal knitr function.
#
#' @keywords internal
escapeLatex = function(x, newlines = FALSE, spaces = FALSE) {
  x = gsub('\\\\', '\\\\textbackslash', x)
  x = gsub('([#$%&_{}])', '\\\\\\1', x)
  x = gsub('\\\\textbackslash', '\\\\textbackslash{}', x)
  x = gsub('~', '\\\\textasciitilde{}', x)
  x = gsub('\\^', '\\\\textasciicircum{}', x)
  if (newlines) x = gsub('(?<!\n)\n(?!\n)', '\\\\\\\\', x, perl = TRUE)
  if (spaces) x = gsub('  ', '\\\\ \\\\ ', x)
  x
}

# Vectors for renaming columns in varcomp and effdim tables.
renameFrom <- c("genotype", "repId", "rowId", "colId", "subBlock",
                "repId:rowId", "repId:colId", "repId:subBlock",
                "colCoord", "rowCoord", "rowCoordcolCoord",
                "f(colCoord)", "f(rowCoord)",
                "f(colCoord):rowCoord",
                "colCoord:f(rowCoord)",
                "f(colCoord):f(rowCoord)", "Nobs", "R", "variance",
                "pow", "units")
renameTo <- c("Genotype", "Replicate", "Row", "Col", "Block",
              "Row(replicate)", "Col(replicate)", "Block(replicate)",
              "Linear trend along cols", "Linear trend along rows",
              "Linear trend along rows and cols",
              "Smooth trend along cols", "Smooth trend along rows",
              "Linear trend in rows changing smoothly along cols",
              "Linear trend in cols changing smoothly along rows",
              "Smooth-by-smooth interaction trend over rows and cols",
              "Number of observations", "Residual", "Residual", "Power",
              "Units")

#' Function for extracting the table with variance components from a model in
#' a nicely printable format.
#'
#' @keywords internal
extractVarComp <- function(model,
                           engine) {
  if (engine == "SpATS") {
    ## Extract variance components directly from model since using summary
    ## creates a matrix with values already rounded restricting flexibility.
    varComp <- matrix(data = c(model$var.comp, model$psi[1]),
                      dimnames = list(c(names(model$var.comp), "Residual"),
                                      "Variance"))

  } else if (engine == "lme4") {
    if (inherits(model, "lm")) {
      ## In this case there is only residual variance since there are no
      ## random effects.
      varComp <- matrix(data = summary(model)$sigma ^ 2,
                        dimnames = list("Residual", "Variance"))
    } else {
      varComp <- as.data.frame(lme4::VarCorr(model))
      varComp <- matrix(data = varComp$vcov,
                        dimnames = list(varComp$grp, "Variance"))
    }
  } else if (engine == "asreml") {
    ## asreml provides the SE of the variance components as standard output.
    ## This is included in varComp.
    varComp <- as.matrix(summary(model)$varcomp[c("component", "std.error")])
    ## Remove correlations from output. These are present for spatials models.
    varComp <- varComp[!grepl(pattern = ".cor", x = rownames(varComp)), ,
                       drop = FALSE]
    ## To extract user readable names similar to the other engines from
    ## asreml output split the rownames on "." and "!" Depending on the first
    ## part of the split use the appropriate part as row name.
    rownames(varComp) <- sapply(X = strsplit(x = rownames(varComp),
                                             split = "[!.]+"),
                                FUN = function(split) {
                                  if (split[[1]] == "R") {
                                    return(split[[2]])
                                  } else {
                                    return(split[[1]])
                                  }
                                })
    colnames(varComp) <- c("Variance", "SE")
  }
  ## Rename rows for more user readable output.
  for (j in seq_along(renameFrom)) {
    rownames(varComp)[rownames(varComp) == renameFrom[j]] <- renameTo[j]
  }
  ## Always put genotype as first row.
  if ("Genotype" %in% rownames(varComp)) {
    varComp <- rbind(varComp["Genotype", , drop = FALSE],
                     varComp[rownames(varComp) != "Genotype", , drop = FALSE])
  }
  ## Add an empty row before residuals if it is not already there.
  ## Only done if there is more than 1 row.
  resRow <- which(rownames(varComp) == "Residual")
  if (nrow(varComp) > 1 && rownames(varComp)[resRow - 1] != "") {
    varComp <- rbind(varComp[1:(resRow - 1), , drop = FALSE],
                     rep(NA, times = ncol(varComp)),
                     varComp[resRow:nrow(varComp), , drop = FALSE])
  }
  return(varComp)
}

#' Helper function for printing anova table in reports.
#' @keywords internal
printAnova <- function(aovTab,
                       title = NULL) {
  ## Add significance stars
  aovTab[, ncol(aovTab) + 1] <-
    symnum(x = aovTab[, ncol(aovTab)], corr = FALSE, na = FALSE,
           cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
           symbols = c("***", "**", "*", ".", " "))
  colnames(aovTab)[ncol(aovTab)] <- ""
  legendText <- paste("Significance codes:",
                      attr(x = aovTab[, ncol(aovTab)], which = "legend"),
                      "} \\\\")
  print(xtable::xtable(x = aovTab, caption = title,
                       label = paste0("anova", title),
                       align = c("l", "r", "r", "r", "r", "r", "l"),
                       digits = c(0, 0, 0, 0, 2, -2, 0),
                       display = c("s", "f", "f", "f", "f", "e", "s")),
        caption.placement = "top",
        latex.environments = "flushleft",
        include.rownames = TRUE, include.colnames = TRUE,
        add.to.row = list(pos = list(nrow(aovTab)),
                          command = paste0("\\hline  \\multicolumn{",
                                           ncol(aovTab), "}{c}{", legendText)))
}


calcPlotBorders <- function(trDat,
                            bordVar) {
  yMin <- min(trDat$rowCoord)
  yMax <- max(trDat$rowCoord)
  xMin <- min(trDat$colCoord)
  xMax <- max(trDat$colCoord)
  ## Create matrix containing replicates.
  ## First create an empty matrix containing all row/column values
  ## between min and max to assure complete missing rows/columns
  ## are added.
  M <- matrix(nrow = yMax - yMin + 1, ncol = xMax - xMin + 1,
              dimnames = list(yMin:yMax, xMin:xMax))
  for (i in 1:nrow(trDat)) {
    M[as.character(trDat[i, "rowCoord"]),
      as.character(trDat[i, "colCoord"])] <- trDat[i, bordVar]
  }
  ## Create an imputed version of M for plotting borders around NA values.
  MImp <- M
  MImp[is.na(MImp)] <- nlevels(trDat[[bordVar]]) + 1
  has.breaks <- function(x) {
    ncol(x) == 2 & nrow(x) > 0
  }
  ## Create a data.frame with positions where the value of rep in the
  ## data changes in vertical direction.
  vertW <- do.call(rbind.data.frame,
                   Filter(f = has.breaks, x = Map(function(i, x) {
                     cbind(y = i, x = which(diff(c(0, x, 0)) != 0))
                   }, 1:nrow(MImp), split(MImp, 1:nrow(MImp)))))
  ## Remove vertical walls that are on the outside bordering an NA value
  ## to prevent drawing of unneeded lines.
  vertW <- vertW[!(vertW$x == 1 & is.na(M[vertW$y, 1])) &
                   !(vertW$x == ncol(M) + 1 &
                       is.na(M[vertW$y, ncol(M)])), ]
  ## Add min row value for plotting in the correct position.
  vertW$y <- vertW$y + yMin - 1
  vertW$x <- vertW$x + xMin - 1
  ## For horizontal walls follow the same procedure as above.
  horW <- do.call(rbind.data.frame,
                  Filter(f = has.breaks, x = Map(function(i, y) {
                    cbind(x = i, y = which(diff(c(0, y, 0)) != 0))
                  }, 1:ncol(MImp), as.data.frame(MImp))))
  horW <- horW[!(horW$y == 1 & is.na(M[1, horW$x])) &
                 !(horW$y == nrow(M) + 1 & is.na(M[nrow(M), horW$x])), ]
  horW$y <- horW$y + yMin - 1
  horW$x <- horW$x + xMin - 1
  return(list(horW = horW, vertW = vertW))
}

## This function is a slightly modified copy of map_data from ggplot2 combined
## with map.fortify also from ggplot2.
## Using the normal function is not possible because both qtl and maps have
## a class map and when building the vignette this gives an error.
mapData <- function(xLim,
                    yLim) {
  mapObj <- maps::map("world", exact = FALSE, plot = FALSE,
                       fill = TRUE, xlim = xLim, ylim = yLim, projection = )
  df <- data.frame(long = mapObj$x, lat = mapObj$y)
  df$group <- cumsum(is.na(df$long) & is.na(df$lat)) + 1
  df$order <- 1:nrow(df)
  names <- do.call("rbind", lapply(strsplit(mapObj$names, "[:,]"),
                                   "[", 1:2))
  df$region <- names[df$group, 1]
  df$subregion <- names[df$group, 2]
  return(df[stats::complete.cases(df$lat, df$long), ])
}
