#' Fits a mixed model.
#'
#' This function fits a mixed model best fitting to the data in a TD object.
#' The exact model fitted is determined by both the data and the chosen
#' parameters.\cr\cr
#' Six different types of models can be fitted depending on the environmental
#' dimension of the data. The enviromental demensions are described below
#' together with their corresponding models and which parameters to specify to
#' fit them.
#' \itemize{
#' \item{environments correspond to trials\cr
#' trait = trial + \strong{genotype + genotype:trial}\cr
#' no extra parameters}
#' \item{trials form a factorial structure of locations x years\cr
#' trait = year + location + year:location + \strong{genotype + genotype:year +
#' genotype:location + genotype:year:location}\cr
#' set \code{locationYear = TRUE}}
#' \item{trials are nested within year\cr
#' trait = year + year:trial + \strong{genotype + genotype:year +
#' genotype:year:trial}\cr
#' set \code{nesting = "year"}}
#' \item{trials are nested within locations\cr
#' trait = location + location:trial + \strong{genotype + genotype:location +
#' genotype:location:trial}\cr
#' set \code{nesting = "loc"}}
#' \item{trials correspond to locations within regions across years\cr
#' trait = region + region:location + year + region:year + region:location:year +
#' \strong{genotype + genotype:region:location + genotype:year +
#' genotype:region:year + genotype:region:location:year}\cr
#' set \code{regionLocatoinYear = TRUE}}
#' \item{trials are nested within scenarios\cr
#' trait = scenario + scenario:trial + \strong{genotype + genotype:scenario +
#' genotype:scenario:trial}\cr
#' set \code{nesting = "scenario"}}
#' }
#' In the models above the random part of the model is printed bold.\cr
#' Note that the final random term is only fitted if the data contains
#' replicates or weights are applied (using the option \code{useWt}).\cr\cr
#' The function first fits a model where all model terms are included as fixed
#' terms. Based on the anova table of this model, terms in the fixed part of the
#' model that are likely to give a problem when fitting the mixed model are
#' removed. Also a warning is printed if the mean sum of squares for a model
#' term points to a possible zero variance component in the mixed model.\cr\cr
#' Then a model is fitted where all model terms are included as random terms.
#' Based on the variance components in this model the percentage of variance
#' explainded by each of the model components is determined. This is printed in
#' the model summary.\cr\cr
#' Finally a mixed model is fitted as specified in the overview above. Based on
#' this mixed model variance components can be computed using \code{\link{vc}},
#' heritabilies can be computed using \code{\link{herit}} and predictions can be
#' made using \code{\link{predict.varComp}}.
#'
#' @inheritParams gxeAmmi
#'
#' @param engine A character string specifying the engine used for modeling.
#' Either "lme4" or "asreml".
#' @param locationYear Should a model be fitted assuming a factorial structure
#' of locations x years?
#' @param nesting A character string specifying a column in TD specifying the
#' nesting structure of the trials.
#' @param regionLocationYear Should a model be fitted assuming locations within
#' regions across years?
#' @param useWt Should the model be fitted using weights? Doing so requieres a
#' column wt in the data.
#' @param diagnostics Should diagnostics on missing combinations of model
#' variables be printed?
#'
#' @return An object of class \code{varComp}, a list containing:
#' \item{fitMod}{The fitted model.}
#' \item{modDat}{A data.frame containing the data used when fitting the model.}
#' \item{nesting}{A name of the variable used as nesting variable in the model.}
#' \item{useLocYear}{A boolean specifying if a model containing location x year
#' interaction was fitted.}
#' \item{fullRandVC}{A data.frame containing the variance components for the
#' fully random model.}
#' \item{aovFullMixedMod}{A data.frame containing the anova table for the fully
#' fixed model.}
#' \item{engine}{The engine used for fitting the model.}
#' \item{diagTabs}{A list of data.frame, one for each random model term,
#' containing the missing combinations in the data for that term.}
#'
#' @examples
#' ## Fit a mixed model.
#' geVarComp <- gxeVarComp(TD = TDMaize, trait = "yld")
#'
#' ## Summarize results.
#' summary(geVarComp)
#'
#' ## Plot the standard deviations.
#' plot(geVarComp)
#'
#' @importFrom utils tail
#' @export
gxeVarComp <- function(TD,
                       trials = names(TD),
                       trait,
                       engine = c("lme4", "asreml"),
                       locationYear = FALSE,
                       nesting = NULL,
                       regionLocationYear = FALSE,
                       useWt = FALSE,
                       diagnostics = FALSE) {
  ## Checks.
  if (missing(TD) || !inherits(TD, "TD")) {
    stop("TD should be a valid object of class TD.\n")
  }
  trials <- chkTrials(trials, TD)
  TDTot <- Reduce(f = rbind, x = TD[trials])
  chkCol(trait, TDTot)
  chkCol("trial", TDTot)
  chkCol("genotype", TDTot)
  if (locationYear) {
    chkCol("loc", TDTot)
    chkCol("year", TDTot)
    ## Check whether year x locations are crossed or nested.
    locYearTab <- table(TDTot[["year"]], TDTot[["loc"]])
    locYearCrossed <- all(locYearTab > 0)
  }
  if (!is.null(nesting)) {
    chkCol(nesting, TDTot)
  }
  if (regionLocationYear) {
    chkCol("loc", TDTot)
    chkCol("year", TDTot)
    chkCol("region", TDTot)
  }
  if (useWt) {
    chkCol("wt", TDTot)
  } else {
    TDTot[["wt"]] <- 1
  }
  engine <- match.arg(engine)
  TDTot <- droplevels(TDTot)
  ## Increase maximum number of iterations for asreml.
  maxIter <- 200
  ## Asreml can't handle missing weights, so set them to 0 when missing.
  TDTot[is.na(TDTot[["wt"]]), "wt"] <- 0
  ## Construct formula for fixed part - first as text.
  ## Trying to fit this in something 'smart' actually makes it unreadable.
  ## First create a vector with the separate terms.
  ## This avoids difficult constructions to get the +-es correct.
  fixedTerms <- c(if (!locationYear && is.null(nesting) &&
                      !regionLocationYear) "trial",
                  if (locationYear) c("year", if (locYearCrossed) "loc",
                                      "year:loc"),
                  if (!is.null(nesting)) c(nesting, paste0(nesting, ":trial")),
                  if (regionLocationYear) c("region", "region:loc", "year",
                                            "region:year", "region:loc",
                                            "region:loc:year"))
  ## Check if the data contains replicates.
  repTab <- table(TDTot[c("genotype",
                          unlist(strsplit(x = tail(fixedTerms, 1), split = ":")))])
  hasReps <- any(repTab > 1)
  ## Random terms are genotype x fixedTerms.
  ## If there are no replicates or weights the final random term is the actual
  ## residual and therefore left out of the model.
  if (hasReps || useWt) {
    randTermIncl <- fixedTerms
  } else {
    randTermIncl <- fixedTerms[-length(fixedTerms)]
  }
  randTerms <- c("genotype",
                 if (length(randTermIncl) > 0) paste0("genotype:", randTermIncl))
  ## First fit a model with all terms fixed to determine:
  ## - should all terms in the fixed part really be present.
  ## - Predict which terms in the random part of the model will probably
  ##   have a zero variance component.
  fullFixedTxt <- paste0("`", trait, "`~",
                         paste(c(fixedTerms, randTerms), collapse = "+"))
  ## Fit the fully fixed model.
  fullFixedMod <- suppressWarnings(lm(formula(fullFixedTxt), data = TDTot))
  aovFullFixedMod <- anova(fullFixedMod)
  rownames(aovFullFixedMod)[nrow(aovFullFixedMod)] <- "residuals"
  ## Get all model terms as used by lm (might involve reordered terms).
  fullFixedLabs <- attr(x = terms(fullFixedMod), which = "term.labels")
  if (!all(fullFixedLabs %in% rownames(aovFullFixedMod))) {
    ## At least one terms missing from anova.
    ## If this is a fixed term remove it from fixed.
    missTerms <- fullFixedLabs[!fullFixedLabs %in% rownames(aovFullFixedMod)]
    for (missTerm in missTerms) {
      fixedTermSets <- strsplit(x = fixedTerms, split = ":")
      missTermSet <- unlist(strsplit(x = missTerm, split = ":"))
      remPos <- !sapply(X = fixedTermSets, FUN = setequal, missTermSet)
      fixedTerms <- fixedTerms[remPos]
    }
  }
  ## Get rand terms that indicate zero variance components.
  ## lm reorders variables in model terms.
  ## use sets of variables in terms to compare them.
  aovTermSets <- strsplit(x = rownames(aovFullFixedMod), split = ":")
  for (randTerm in randTerms) {
    ## Convert term to set of variables in term.
    randTermSet <- unlist(strsplit(x = randTerm, split = ":"))
    ## Get position of term in anova table by comparing sets.
    randTermPos <- sapply(X = aovTermSets, FUN = setequal, randTermSet)
    ## Get MSS for current term.
    MSSRandTerm <- aovFullFixedMod[randTermPos, "Mean Sq"]
    ## For all other terms in the anova table that have the current term
    ## as a subset the MSS cannot be higher.
    ## If it is the corresponding variance component is possibly zero.
    for (i in seq_along(aovTermSets)) {
      if ((all(randTermSet %in% aovTermSets[i]) ||
           ## Always include the residual term for comparison.
           i == nrow(aovFullFixedMod)) &&
          !is.nan(aovFullFixedMod[i, "Mean Sq"]) &&
          aovFullFixedMod[i, "Mean Sq"] > MSSRandTerm) {
        warning("Mean Sum of Squares for ", randTerm, " smaller than Mean ",
                "Sum of Squares for ", rownames(aovFullFixedMod)[i], ".\n",
                "Possible zero variance components.\n", call. = FALSE)
      }
    }
  }
  ## Fit the fully random model and compute how much variation
  ## is explained by each of the model terms.
  ## This is stored as fullRandVC and included in the output to create
  ## a nice summary.
  if (engine == "lme4") {
    ## Construct input for full random model.
    fullRandTxt <- paste0("`", trait, "`~",
                          paste(paste0("(1|", c(fixedTerms, randTerms), ")"),
                                collapse = "+"))
    fullRandMod <- lme4::lmer(formula(fullRandTxt), data = TDTot)
    fullRandVC <- as.data.frame(lme4::VarCorr(fullRandMod))
    rownames(fullRandVC) <- fullRandVC[["grp"]]
    rownames(fullRandVC)[nrow(fullRandVC)] <- "residuals"
    vcovTot <- sum(fullRandVC[["vcov"]])
    fullRandVC[["vcovPerc"]] <- fullRandVC[["vcov"]] / vcovTot
    fullRandVC <- fullRandVC[c((nrow(fullRandVC)-1):1, nrow(fullRandVC)),
                             c("vcov", "vcovPerc")]
  } else if (engine == "asreml") {
    ## Construct input for full random model.
    fullRandTxt <- paste("~", paste(c(fixedTerms, randTerms), collapse = "+"))
    fullRandMod <- tryCatchExt(asreml::asreml(fixed = formula(paste0("`", trait, "`~ 1")),
                                              random = formula(fullRandTxt), data = TDTot,
                                              trace = FALSE))
    if (!is.null(fullRandMod$warning)) {
      ## Check if param 1% increase is significant. Remove warning if not.
      fullRandMod <- chkLastIter(fullRandMod)
    }
    if (length(fullRandMod$warning) != 0) {
      warning(fullRandMod$warning, "\n", call. = FALSE)
    }
    if (!is.null(fullRandMod$error)) {
      stop(fullRandMod$error)
    } else {
      fullRandMod <- fullRandMod$value
    }
    fullRandVC <- summary(fullRandMod)$varcomp
    rownames(fullRandVC)[nrow(fullRandVC)] <- "residuals"
    vcovTot <- sum(fullRandVC[["component"]])
    fullRandVC[["vcovPerc"]] <- fullRandVC[["component"]] / vcovTot
    colnames(fullRandVC)[colnames(fullRandVC) == "component"] <- "vcov"
    fullRandVC <- fullRandVC[, c("vcov", "vcovPerc")]
  }
  ## Create tables for diagnostics.
  diagTabs <- lapply(X = fixedTerms, FUN = function(fixedTerm) {
    fixedVars <- unlist(strsplit(x = fixedTerm, split = ":"))
    fixedVarLevs <- lapply(X = fixedVars, FUN = function(fixedVar) {
      unique(TDTot[[fixedVar]])
    })
    fullTab <- expand.grid(c(list(unique(TDTot[["genotype"]])), fixedVarLevs),
                           KEEP.OUT.ATTRS = FALSE)
    missTab <- fullTab[!interaction(fullTab) %in%
                         interaction(TDTot[c("genotype", fixedVars)]), ]
    colnames(missTab) <- c("genotype", fixedVars)
    return(missTab)
  })
  ## Create the full fixed part of the model as a character.
  ## This is identical for lme4 and asreml so only needs to be done once.
  fixedTxt <- paste0("`", trait, "`~", paste(fixedTerms, collapse = "+"))
  if (engine == "lme4") {
    randTxt <- paste(paste0("(1|", randTerms, ")"), collapse = "+")
    formTxt <- paste(fixedTxt, "+", randTxt)
    ## Fit the actual model.
    mr <- lme4::lmer(formula(formTxt), data = TDTot, weights = TDTot[["wt"]])
    ## Construct STA object.
  } else if (engine == "asreml") {
    if (requireNamespace("asreml", quietly = TRUE)) {
      randTxt <- paste("~ ", paste(randTerms, collapse = "+"))
      ## Put arguments for models in a list to make it easier to switch
      ## between asreml3 and asreml4. Usually only one or two arguments differ.
      ## Also some arguments are identical for all models
      modArgs0 <- list(fixed = formula(fixedTxt), random = formula(randTxt),
                       data = TDTot, weights = "wt", maxiter = maxIter,
                       trace = FALSE)
      modArgs <- modArgs0
      ## Fit the actual model.
      mr <- tryCatchExt(do.call(asreml::asreml, modArgs))
      if (!is.null(mr$warning)) {
        ## Check if param 1% increase is significant. Remove warning if not.
        mr <- chkLastIter(mr)
      }
      if (length(mr$warning) != 0) {
        warning("Asreml gave the following warning:\n", mr$warning, "\n",
                call. = FALSE)
      }
      if (!is.null(mr$error)) {
        warning("Asreml gave the following error:\n", mr$error, call. = FALSE)
      } else {
        mr <- mr$value
        mr$call$fixed <- eval(mr$call$fixed)
        mr$call$random <- eval(mr$call$random)
      }
    } else {
      stop("Failed to load 'asreml'.\n")
    }
  }
  if (diagnostics) {
    ## Print diagnostics.
    for (diagTab in diagTabs) {
      if (nrow(diagTab) > 0) {
        cat(nrow(diagTab), " missing combinations for ",
            paste(colnames(diagTab), collapse = " x "), ".\n", sep = "")
        if (nrow(diagTab) > 10) {
          cat("Printing first 10 missing combinations.\n",
              "Use diagnostics() on the result for the full data.frame.\n")
        }
        print(head(diagTab, 10), row.names = FALSE)
        cat("\n\n")
      } else {
        cat("No missing combinations for ",
            paste(colnames(diagTab), collapse = " x "), ".\n\n", sep = "")
      }
    }
  }
  ## Create output.
  res <- createVarComp(fitMod = mr, modDat = TDTot, nesting = nesting,
                       useLocYear = locationYear, fullRandVC = fullRandVC,
                       aovFullFixedMod = aovFullFixedMod, engine = engine,
                       diagTabs = diagTabs)
  return(res)
}


