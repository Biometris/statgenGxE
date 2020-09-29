#' Helper function for suppressing a single warning message.
#'
#' @noRd
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
#' If the corresponding parameter is close to zero and then changes of 1%
#' or more can be expected and are ok.
#'
#' @noRd
#' @keywords internal
chkLastIter <- function(model) {
  wrnMsg <- "changed by more than 1%"
  if (any(grepl(pattern = wrnMsg, x = model$warning))) {
    if (asreml4()) {
      ## EXtract trace df from model object.
      mon <- model$value$trace
      ## Extract values for parameters for last 2 iterations.
      ## First 3 rows give general model info.
      lastIt <- mon[-(1:3), c(ncol(mon) - 1, ncol(mon))]
    } else {
      ## EXtract monitor df from model object.
      mon <- model$value$monitor
      ## Extract values for parameters for last 2 iterations.
      ## First 3 rows give general model info. Last col a summary.
      lastIt <- mon[-(1:3), c(ncol(mon) - 2, ncol(mon) - 1)]
    }
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
#'
#' @noRd
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
#' @noRd
#' @keywords internal
predictAsreml <- function(model,
                          classify = "genotype",
                          associate = as.formula("~ NULL"),
                          vcov = TRUE,
                          TD,
                          ...) {
  wrnMsg <- "reset workspace or pworkspace arguments"
  ## Predict using default settings, i.e. pworkspace = 8e6
  modelP <- tryCatchExt(predict(model, classify = classify,
                                vcov = vcov, associate = associate,
                                data = TD, maxiter = 20, trace = FALSE, ...))
  pWorkSpace <- 8e6
  ## While there is a warning, increase pWorkSpace and predict again.
  while (!is.null(modelP$warning) &&
         any(grepl(pattern = wrnMsg, x = modelP$warning))
         && pWorkSpace < 160e6) {
    pWorkSpace <- pWorkSpace + 8e6
    modelP <- tryCatchExt(predict(model, classify = classify,
                                  vcov = vcov, associate = associate, data = TD,
                                  maxiter = 20, pworkspace = pWorkSpace,
                                  trace = FALSE, ...))
  }
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
    stop("Error in asreml when running predict. Asreml message:\n",
         modelP$error, "\n",
         modelP$warning, "\n", call. = FALSE)
  }
}

#' Helper function for printing anova table in reports.
#'
#' @noRd
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

#' Function for escaping special LaTeX characters
#'
#' Taken from knitr package. Copied since it is an internal knitr function.
#'
#' @noRd
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

#' Helper function for detecting the version of asreml installed.
#' This is used wherever the syntax for asreml4 differs from asreml3.
#'
#' @noRd
#' @importFrom utils packageVersion
#' @keywords internal
asreml4 <- function() {
  if (requireNamespace("asreml", quietly = TRUE)) {
    if (packageVersion("asreml")[1] >= 4) {
      ## Calling license status apparently also activates the license if this
      ## was done once before.
      licenceStatus <- asreml::asreml.license.status(quiet = TRUE)
      if (licenceStatus$status != 0) {
        stop("Error checking asreml licence status:\n",
             licenceStatus$statusMessage)
      }
      return(TRUE)
    }
    return(FALSE)
  }
}
