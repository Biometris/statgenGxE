#' S3 class AMMI
#'
#' Function for creating objects of S3 class AMMI.\cr
#' \code{\link{print}}, \code{\link{summary}}, \code{\link{plot}} and
#' \code{\link{report}} methods are available.
#'
#' @param envScores A matrix containing environmental scores.
#' @param genoScores A matrix containing genotypic scores.
#' @param importance A data.frame containing the importance of the principal
#' components.
#' @param anova A data.frame containing anova scores of the AMMI analysis.
#' @param fitted A matrix containing fitted values from the AMMI model.
#' @param trait A character string indicating the analyzed trait.
#' @param envMean A numerical vector containing the environmental means.
#' @param genoMean A numerical vector containing the genotypic means.
#' @param overallMean A numerical value containing the overall mean.
#' @param GGE Has a GGE analysis been performed?
#' @param byYear Has the analysis been performed by year?
#'
#' @seealso \code{\link{plot.AMMI}}, \code{\link{report.AMMI}}
#'
#' @name AMMI
NULL

#' @rdname AMMI
#' @keywords internal
createAMMI <- function(envScores,
                       genoScores,
                       importance,
                       anova,
                       fitted,
                       trait,
                       envMean,
                       genoMean,
                       overallMean,
                       dat,
                       GGE,
                       byYear) {
  AMMI <- structure(list(envScores = envScores,
                         genoScores = genoScores,
                         importance = importance,
                         anova = anova,
                         fitted = fitted,
                         trait = trait,
                         envMean = envMean,
                         genoMean = genoMean,
                         overallMean = overallMean,
                         dat = dat,
                         GGE = GGE,
                         byYear = byYear),
                    class = "AMMI",
                    timestamp = Sys.time())
  return(AMMI)
}

#' @export
print.AMMI <- function(x,
                       ...,
                       printGenoScores = FALSE) {
  cat("Principal components",
      "\n====================\n")
  if (x$byYear) {
    years <- names(x$importance)
    for (year in years) {
      cat(paste(year, "\n"))
      print(x$importance[[year]][, 1:ncol(x$envScores[[year]]), drop = FALSE])
      cat("\n")
    }
  } else {
    print(x$importance[, 1:ncol(x$envScores), drop = FALSE])
  }
  cat("\nAnova",
      "\n=====\n")
  if (x$byYear) {
    for (year in years) {
      cat(paste(year, "\n"))
      print(x$anova[[year]])
      cat("\n")
    }
  } else {
    print(x$anova)
  }
  cat("\nEnvironment scores",
      "\n==================\n")
  if (x$byYear) {
    for (year in years) {
      cat(paste(year, "\n"))
      print(x$envScores[[year]], ...)
      cat("\n")
    }
  } else {
    print(x$envScores, ...)
  }
  if (printGenoScores) {
    cat("\nGenotypic scores",
        "\n================\n")
    if (x$byYear) {
      for (year in years) {
        cat(paste(year, "\n"))
        print(x$genoScores[[year]], ...)
        cat("\n")
      }
    } else {
      print(x$genoScores, ...)
    }
  }
}

#' @export
summary.AMMI <- function(object,
                         ...,
                         printGenoScores = FALSE) {
  print(object, ..., printGenoScores = printGenoScores)
}

#' Plot function for class AMMI
#'
#' Two types of biplot can be made. A plot of genotype and environment
#' means vs PC1 (AMMI1) or a biplot of genotypes and environment interaction
#' with PC1 and PC2 (AMMI2).\cr\cr
#' If the AMMI analysis was done by year, a separate plot will be made for
#' every year in the data. For some years the number of principal components
#' may be lower than the number specified on the secondary axis. If this is the
#' case this year is skipped when plotting. If this happens for all years the
#' function returns an error.
#'
#' @param x An object of class AMMI
#' @param ... Not used.
#' @param plotType A character string indicating which plot should be made.
#' Either "AMMI1" for an AMMI1 plot (genotype and
#' environment means vs PC1) or "AMMI2" for an AMMI2 biplot (genotypes and
#' environment interaction with PC1 and PC2) respectively. For results of a
#' GGE analysis only an GGE2 biplot can be made and plotType may be ignored.
#' @param scale A numerical value. The variables are scaled by
#' \code{lambda ^ scale} and the observations by \code{lambda ^ (1 - scale)}
#' where \code{lambda} are the singular values computed by
#' \code{\link[stats]{princomp}} in \code{\link{gxeAmmi}}. Normally
#' \code{0 <= scale <= 1}, and a warning will be issued if the specified
#' scale is outside this range.
#' @param plotGeno Should genotypes be plotted?
#' @param colorGenoBy A character string indicating a column in the \code{TD}
#' used as input for the AMMI analysis by which the genotypes should be colored.
#' If \code{NULL} all genotypes will be colored in black.
#' @param colGeno A character vector with plot colors for the genotypes. A
#' single color when \code{colorGenoBy = NULL}, a vector of colors otherwise.
#' @param sizeGeno An numerical value indicating the text size for plotting the
#' genotypes. Use \code{sizeGeno = 0} for plotting genotypes as points instead
#' of using their names.
#' @param plotConvHull Should a convex hull be plotted around the genotypes. If
#' \code{TRUE} a convex hull is plotted. For GGE2 biplots lines from the origin
#' of the plot perpendicular to the edges of the hull are added. Only valid for
#' AMMI2 and GGE2 biplots.
#' @param plotEnv Should environments be plotted?
#' @param colorEnvBy A character string indicating a column in the \code{TD}
#' used as input for the AMMI analysis by which the environments should be
#' colored. If \code{NULL} all genotypes will be colored in red.
#' @param colEnv A character string with plot colors for the environments. A
#' single color when \code{colorEnvBy = NULL}, a vector of colors otherwise.
#' @param sizeEnv An integer indicating the text size for plotting the
#' environments.
#' @param envFactor A positive numerical value giving a factor by which to blow
#' up the environmental scores. Providing a value between 0 and 1 will
#' effectively blow up the genotypic scores.
#' @param primAxis A character string indicating the principal component to be
#' plotted on the primary axis of the AMMI2 plot. Has to be given as
#' \code{"PCn"} where n is the number of the principal component.
#' @param secAxis A character string indicating the principal component to be
#' plotted on the secondary axis of the AMMI2 plot. Has to be given as
#' \code{"PCn"} where n is the number of the principal component. n Has to
#' differ from \code{primAxis}.
#' @param rotatePC A character string indicating a genotype or environment that
#' is to be aligned with the positive x-axis in the plot.
#' @param title A character string used a title for the plot.
#' @param output Should the plot be output to the current device? If
#' \code{FALSE} only a list of ggplot objects is invisibly returned.
#'
#' @return A biplot depending on \code{plotType}. The ggplot object for the
#' biplot is returned invisibly.
#'
#' @examples
#' ## Run AMMI analysis.
#' geAmmi <- gxeAmmi(TD = TDMaize, trait = "yld")
#'
#' ## Create an AMMI1 biplot.
#' plot(geAmmi)
#'
#' ## Create an AMMI2 biplot.
#' plot(geAmmi, plotType = "AMMI2", scale = 0.5)
#'
#' ## Create an AMMI2 biplot, with HN96b along the positive x-axis.
#' plot(geAmmi, plotType = "AMMI2", scale = 0.5, rotatePC = "HN96b")
#'
#' ## Run GGE analysis.
#' geGGE <- gxeGGE(TD = TDMaize, trait = "yld")
#'
#' ## Create an GGE2 biplot.
#' ## Add a convex hull.
#' plot(geGGE, plotType = "GGE2", scale = 0.5, plotConvHull = TRUE)
#'
#' @importFrom grDevices topo.colors
#'
#' @family AMMI
#'
#' @export
plot.AMMI <- function(x,
                      ...,
                      plotType = c("AMMI1", "AMMI2", "GGE2"),
                      scale = 1,
                      plotGeno = TRUE,
                      colorGenoBy = NULL,
                      colGeno = NULL,
                      sizeGeno = 0,
                      plotConvHull = FALSE,
                      plotEnv = TRUE,
                      colorEnvBy = NULL,
                      colEnv = NULL,
                      sizeEnv = 3,
                      envFactor = 1,
                      primAxis = "PC1",
                      secAxis = "PC2",
                      rotatePC = NULL,
                      title = NULL,
                      output = TRUE) {
  ## Checks.
  if (missing(plotType) && x$GGE) {
    plotType <- "AMMI2"
  } else {
    plotType <- match.arg(plotType)
    if (plotType == "GGE2") {
      plotType <- "AMMI2"
    }
  }
  chkChar(title, len = 1)
  chkNum(scale, min = 0, max = 1, null = FALSE, incl = TRUE)
  if (plotGeno) {
    chkNum(sizeGeno, min = 0, null = FALSE, incl = TRUE)
    chkChar(colGeno)
  }
  if (plotEnv) {
    chkNum(sizeEnv, min = 0, null = FALSE, incl = TRUE)
    chkChar(colEnv)
  }
  chkChar(colorGenoBy)
  if (!is.null(colorGenoBy)) {
    if (x$byYear) {
      for (dat in x$dat) {
        chkCol(colorGenoBy, dat)
        colTab <- unique(dat[c("genotype", colorGenoBy)])
        if (nrow(colTab) != nlevels(droplevels(dat[["genotype"]]))) {
          stop("colorGenoBy should have exactly one value per genotype.\n")
        }
        if (any(is.na(dat[[colorGenoBy]]))) {
          stop("Missing values in ", colorGenoBy, ".\n")
        }
      }
    } else {
      chkCol(colorGenoBy, x$dat)
      colTab <- unique(x$dat[c("genotype", colorGenoBy)])
      if (any(is.na(x$dat[[colorGenoBy]]))) {
        stop("Missing values in ", colorGenoBy, ".\n")
      }
      if (nrow(colTab) != nlevels(droplevels(x$dat[["genotype"]]))) {
        stop("colorGenoBy should have exactly one value per genotype.\n")
      }
    }
  }
  chkChar(colorEnvBy)
  if (!is.null(colorEnvBy)) {
    if (x$byYear) {
      for (dat in x$dat) {
        chkCol(colorEnvBy, dat)
        colTab <- unique(dat[c("trial", colorEnvBy)])
        if (nrow(colTab) != nlevels(droplevels(dat[["trial"]]))) {
          stop("colorEnvBy should have exactly one value per environment.\n")
        }
        if (any(is.na(dat[[colorEnvBy]]))) {
          stop("Missing values in ", colorEnvBy, ".\n")
        }
      }
    } else {
      chkCol(colorEnvBy, x$dat)
      colTab <- unique(x$dat[c("trial", colorEnvBy)])
      if (nrow(colTab) != nlevels(droplevels(x$dat[["trial"]]))) {
        stop("colorEnvBy should have exactly one value per environment.\n")
      }
      if (any(is.na(x$dat[[colorEnvBy]]))) {
        stop("Missing values in ", colorEnvBy, ".\n")
      }
    }
  }
  chkNum(envFactor, min = 0)
  chkChar(rotatePC, len = 1)
  if (!is.null(rotatePC) && !rotatePC %in% x$dat[["genotype"]] &&
      !rotatePC %in% x$dat[["trial"]]) {
    stop("rotatePC should specify a genotype or trial present in data.\n")
  }
  if (plotType == "AMMI1") {
    if (x$byYear) {
      ## Create a list of AMMI1 plots.
      p <- sapply(X = names(x$envScores), FUN = function(year) {
        plotAMMI1(loadings = x$envScores[[year]] * envFactor,
                  scores = x$genoScores[[year]] / envFactor,
                  importance = x$importance[[year]],
                  overallMean = x$overallMean[[year]] * envFactor,
                  genoMean = x$genoMean[[year]] / envFactor,
                  envMean = x$envMean[[year]] * envFactor,
                  trait = x$trait, dat = x$dat[[year]], GGE = x$GGE,
                  year = year, rotatePC = rotatePC, scale = scale,
                  plotGeno = plotGeno, colGeno = colGeno, sizeGeno = sizeGeno,
                  plotEnv = plotEnv, colorEnvBy = colorEnvBy, colEnv = colEnv,
                  sizeEnv = sizeEnv, colorGenoBy = colorGenoBy, title = title)
      }, simplify = FALSE)
    } else {
      ## Create a single AMMI1 plot.
      p <- plotAMMI1(loadings = x$envScores * envFactor,
                     scores = x$genoScores / envFactor,
                     importance = x$importance,
                     overallMean = x$overallMean * envFactor,
                     genoMean = x$genoMean / envFactor,
                     envMean = x$envMean * envFactor,
                     trait = x$trait, dat = x$dat, GGE = x$GGE,
                     rotatePC = rotatePC, scale = scale,
                     plotGeno = plotGeno, colGeno = colGeno,
                     sizeGeno = sizeGeno, plotEnv = plotEnv,
                     colorEnvBy = colorEnvBy,
                     colEnv = colEnv, sizeEnv = sizeEnv,
                     colorGenoBy = colorGenoBy, title = title)
    }
  } else if (plotType == "AMMI2") {
    if (!is.character(primAxis) || length(primAxis) > 1 ||
        substring(text = primAxis, first = 1, last = 2) != "PC") {
      stop("primAxis should be a single character string starting with PC.\n")
    }
    nPC1 <- suppressWarnings(as.numeric(substring(text = primAxis, first = 3)))
    if (is.na(nPC1)) {
      stop("Invalid value provided for primAxis. Make sure the value is ",
           "of the form primAxis = 'PCn' where n is the principal ",
           "component to plot on the primary axis.\n")
    }
    if (!is.character(secAxis) || length(secAxis) > 1 ||
        substring(text = secAxis, first = 1, last = 2) != "PC") {
      stop("secAxis should be a single character string starting with PC.\n")
    }
    nPC2 <- suppressWarnings(as.numeric(substring(text = secAxis, first = 3)))
    if (is.na(nPC2)) {
      stop("Invalid value provided for secAxis. Make sure the value is ",
           "of the form secAxis = 'PCn' where n is the principal ",
           "component to plot on the secondary axis.\n")
    }
    if (nPC1 == nPC2) {
      stop("primAxis should differ from secAxis.\n")
    }
    if (x$byYear) {
      nPCs <- sapply(X = x$envScores, FUN = ncol)
      maxPC <- max(nPCs)
      if (nPC1 > maxPC) {
        stop("Highest number of principal components is ", maxPC,
             ". Plotting of PC", nPC1, " is not possible.\n")
      }
      if (nPC2 > maxPC) {
        stop("Highest number of principal components is ", maxPC,
             ". Plotting of PC", nPC2, " is not possible.\n")
      }
      p <- sapply(X = names(x$envScores)[nPCs >= max(nPC1, nPC2)],
                  FUN = function(year) {
                    plotAMMI2(loadings = x$envScores[[year]] * envFactor,
                              scores = x$genoScores[[year]],
                              importance = x$importance[[year]],
                              trait = x$trait, dat = x$dat[[year]], GGE = x$GGE,
                              year = year, primAxis = primAxis,
                              rotatePC = rotatePC,
                              secAxis = secAxis, scale = scale,
                              plotGeno = plotGeno, colGeno = colGeno,
                              sizeGeno = sizeGeno, plotConvHull = plotConvHull,
                              plotEnv = plotEnv, colEnv = colEnv,
                              colorEnvBy = colorEnvBy, sizeEnv = sizeEnv,
                              colorGenoBy = colorGenoBy, title = title)
                  }, simplify = FALSE)


    } else {
      if (nPC1 > ncol(x$envScores)) {
        stop("AMMI was run with ", ncol(x$envScores), " principal ",
             "components. Plotting of PC", nPC1, " is not possible.\n")
      }
      if (nPC2 > ncol(x$envScores)) {
        stop("AMMI was run with ", ncol(x$envScores), " principal ",
             "components. Plotting of PC", nPC2, " is not possible.\n")
      }
      ## Create a single AMMI2 plot.
      p <- plotAMMI2(loadings = x$envScores * envFactor,
                     scores = x$genoScores,
                     importance = x$importance, trait = x$trait, dat = x$dat,
                     GGE = x$GGE, primAxis = primAxis, secAxis = secAxis,
                     rotatePC = rotatePC, scale = scale, plotGeno = plotGeno,
                     colGeno = colGeno, sizeGeno = sizeGeno,
                     plotConvHull = plotConvHull,
                     plotEnv = plotEnv, colEnv = colEnv,
                     colorEnvBy = colorEnvBy, sizeEnv = sizeEnv,
                     colorGenoBy = colorGenoBy, title = title)
    }
  }
  if (output) {
    if (x$byYear) {
      lapply(X = p, FUN = plot)
    } else {
      plot(p)
    }
  }
  invisible(p)
}

#' Helper function for creating AMMI1 plot
#'
#' @noRd
#' @keywords internal
plotAMMI1 <- function(loadings,
                      scores,
                      importance,
                      overallMean,
                      genoMean,
                      envMean,
                      trait,
                      dat,
                      GGE,
                      year = "",
                      rotatePC = NULL,
                      scale,
                      plotGeno = TRUE,
                      colorGenoBy,
                      colGeno = NULL,
                      sizeGeno,
                      plotEnv = TRUE,
                      colorEnvBy,
                      colEnv = NULL,
                      sizeEnv,
                      title = NULL) {
  percPC1 <- round(importance[2, 1] * 100, 1)
  ## Calculate lambda scale
  lam <- importance[1, 1] ^ scale
  if (is.null(title)) {
    title <- paste(ifelse(GGE, "GGE", "AMMI1"), "plot for", trait, year)
  }
  if (plotGeno) {
    ## Create data.frame for genotypes.
    genoDat <- data.frame(type = "geno", x = genoMean, y = scores[, 1] / lam)
    if (!is.null(colorGenoBy)) {
      ## Merge the colorGenoBy column to the data.
      genoDat <- merge(genoDat, unique(dat[c("genotype", colorGenoBy)]),
                       by.x = "row.names", by.y = "genotype")
      if (!is.factor(genoDat[[colorGenoBy]])) {
        genoDat[[colorGenoBy]] <- as.factor(genoDat[[colorGenoBy]])
      }
      genoDat <- genoDat[order(genoDat[[colorGenoBy]]), ]
      nColGeno <- nlevels(genoDat[[colorGenoBy]])
      if (length(colGeno) == 0) {
        ## Defaults to black for one color for genotypes.
        ## For more than one colors from statgen.genoColors are used.
        ## Fall back to topo.colors if number of colors in option is too small.
        if (nColGeno == 1) {
          colGeno <- "black"
        } else if (length(getOption("statgen.genoColors")) >= nColGeno) {
          colGeno <- getOption("statgen.genoColors")[1:nColGeno]
        } else {
          colGeno <- topo.colors(nColGeno)
        }
      } else {
        nColGenoArg <- length(colGeno)
        if (nColGenoArg != nColGeno) {
          stop("Number of colors provided doesn't match number of genotype groups:",
               "\n", nColGenoArg, " colors provided, ", nColGeno,
               " groups in data.\n")
        }
      }
      ## Add group variable with contents of colorGenoBy.
      genoDat[[".group"]] <- genoDat[[colorGenoBy]]
      ## Set colors as levels of colorGenoBy.
      levels(genoDat[[colorGenoBy]]) <- colGeno
      rownames(genoDat) <- genoDat[["Row.names"]]
      colnames(genoDat)[colnames(genoDat) == colorGenoBy] <- ".color"
    } else {
      genoDat[[".group"]] <- "genoGroup1"
      genoDat[[".color"]] <-
        factor(if (is.null(colGeno)) "black" else colGeno[1])
    }
    ## Restrict columns so rbinding to envDat is possible.
    genoDat <- genoDat[c("type","x", "y", ".group", ".color")]
    ## Add size.
    genoDat[[".size"]] <- sizeGeno
  } else {
    genoDat <- NULL
  }
  if (plotEnv) {
    ## Create data.frame for environments.
    envDat <- data.frame(type = "env", x = envMean, y = loadings[, 1] * lam)
    if (!is.null(colorEnvBy)) {
      envDat <- merge(envDat, unique(dat[c("trial", colorEnvBy)]),
                      by.x = "row.names", by.y = "trial")
      if (!is.factor(envDat[[colorEnvBy]])) {
        envDat[[colorEnvBy]] <- as.factor(envDat[[colorEnvBy]])
      }
      envDat <- envDat[order(envDat[[colorEnvBy]]), ]
      nColEnv <- nlevels(envDat[[colorEnvBy]])
      if (length(colEnv) == 0) {
        ## Defaults to red for one color for trials.
        ## For more than one colors from statgen.trialColors are used.
        ## Fall back to topo.colors if number of colors in option is too small.
        if (nColEnv == 1) {
          colEnv <- "red"
        } else if (length(getOption("statgen.trialColors")) >= nColEnv) {
          colEnv <- getOption("statgen.trialColors")[1:nColEnv]
        } else {
          colEnv <- topo.colors(nColEnv)
        }
      } else {
        nColEnvArg <- length(colEnv)
        if (nColEnvArg != nColEnv) {
          stop("Number of colors provided doesn't match number of environment ",
               "groups:\n", nColEnvArg, " colors provided, ", nColEnv,
               " groups in data.\n")
        }
      }
      ## Add group variable with contents of colorGenoBy.
      envDat[[".group"]] <- envDat[[colorEnvBy]]
      ## Set colors as levels of colorEnvBy
      levels(envDat[[colorEnvBy]]) <- colEnv
      rownames(envDat) <- envDat[["Row.names"]]
      colnames(envDat)[colnames(envDat) == colorEnvBy] <- ".color"
    } else {
      envDat[[".group"]] <- "envGroup1"
      envDat[[".color"]] <- factor(if (is.null(colEnv)) "red" else colEnv[1])
    }
    ## Restrict columns so rbinding to genoDat is possible.
    envDat <- envDat[c("type", "x", "y", ".group", ".color")]
    ## Add size.
    envDat[[".size"]] <- sizeEnv
  } else {
    envDat <- NULL
  }
  ## Create total data,frame for convenience.
  totDat <- rbind(genoDat, envDat)
  ## Rotate in such a way that the selected environment or genotype
  ## is aligned with the positive x-axis.
  if (!is.null(rotatePC)) {
    ## Rotating around x value of overallMean.
    totDat[["x"]] <- totDat[["x"]] + overallMean
    genoDat[["x"]] <- genoDat[["x"]] + overallMean
    envDat[["x"]] <- envDat[["x"]] + overallMean
    ## Get the angle between the selected env/geno and the positive x-axis.
    xRot <- totDat[rotatePC, "x"]
    yRot <- totDat[rotatePC, "y"]
    theta <- atan2(yRot, xRot)
    ## Rotate clockwise over this angle.
    totDat <- rotatePC(dat = totDat, theta = theta, primAxis = "x",
                       secAxis = "y")
    genoDat <- rotatePC(dat = genoDat, theta = theta, primAxis = "x",
                        secAxis = "y")
    envDat <- rotatePC(dat = envDat, theta = theta, primAxis = "x",
                       secAxis = "y")
    ## Rotating around x value of overallMean.
    totDat[["x"]] <- totDat[["x"]] + overallMean
    genoDat[["x"]] <- genoDat[["x"]] + overallMean
    envDat[["x"]] <- envDat[["x"]] + overallMean
  }
  ## Bind together so everything can be plotted in one go.
  ## This has to be done because only one color legend is allowed.
  ## Split between points and text data.
  if (sizeGeno == 0) {
    pointDat <- genoDat
    textDat <- envDat
  } else {
    pointDat <- NULL
    textDat <- rbind(genoDat, envDat)
  }
  ## Create total data,frame for convenience.
  totDat <- rbind(genoDat, envDat)
  ## Get x and y limits and compute plotRatio from them.
  ## This assures 1 unit on the x-axis is identical to 1 unit on the y-axis.
  yMin <- min(totDat[["y"]])
  yMax <- max(totDat[["y"]])
  xMin <- min(totDat[["x"]])
  xMax <- max(totDat[["x"]])
  plotRatio <- (xMax - xMin) / (yMax - yMin)
  colGroups <- setdiff(unique(totDat[[".group"]]),
                       c("envGroup1", "genoGroup1"))
  colNamed <- setNames(c(unique(as.character(genoDat[[".color"]])),
                         unique(as.character(envDat[[".color"]]))),
                       unique(as.character(totDat[[".group"]])))
  if (plotGeno && !is.null(colorGenoBy) && plotEnv && !is.null(colorEnvBy)) {
    nGenoGroups <- length(unique(genoDat[[".group"]]))
    nEnvGroups <- length(unique(envDat[[".group"]]))
    nrowGuide <- nGenoGroups
    shapesGuide <- c(rep(16, times = nGenoGroups), rep(NA, times = nEnvGroups))
    sizesGuide <- c(rep(if (sizeGeno == 0) 2 else sizeGeno,
                        times = nGenoGroups),
                    rep(sizeEnv, times = nEnvGroups))
  } else {
    nrowGuide <- shapesGuide <- sizesGuide <- NULL
  }
  lims <- c(min(xMin, yMin), max(xMax, yMax))
  p <- ggplot2::ggplot(totDat,
                       ggplot2::aes(x = .data[["x"]], y = .data[["y"]],
                                    color = .data[[".group"]])) +
    ## Needed for a square plot output.
    ggplot2::coord_equal(xlim = lims, ylim = lims, clip = "off") +
    ggplot2::scale_color_manual(breaks = colGroups, values = colNamed,
                                name = NULL,
                                guide = ggplot2::guide_legend(nrow = nrowGuide,
                                                              override.aes = list(shape = shapesGuide,
                                                                                  size = sizesGuide))) +
    ## Add reference axes.
    ggplot2::geom_vline(ggplot2::aes(xintercept = overallMean),
                        linetype = "dashed", show.legend = FALSE) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 0),
                        linetype = "dashed", show.legend = FALSE) +
    ## Add labeling.
    ggplot2::labs(x = "Main Effects", y = paste0("PC1 (", percPC1, "%)"),
                  title = title) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   panel.grid = ggplot2::element_blank())
  if (!is.null(pointDat)) {
    ## Plot genotypes as points.
    p <- p + ggplot2::geom_point(data = pointDat,
                                 show.legend = !is.null(colorGenoBy))
  }
  if (!is.null(textDat)) {
    ## Plot genotypes and environments as text.
    p <- p + ggplot2::geom_text(data = textDat,
                                ggplot2::aes(label = rownames(textDat),
                                             size = .data[[".size"]]),
                                show.legend = !is.null(colorEnvBy), vjust = 0) +
      ggplot2::scale_size(range = range(textDat[[".size"]]), guide = "none")
  }
  return(p)
}

#' Helper function for creating AMMI2 plot
#'
#' @noRd
#' @keywords internal
plotAMMI2 <- function(loadings,
                      scores,
                      importance,
                      trait,
                      dat,
                      GGE,
                      year = "",
                      primAxis = "PC1",
                      secAxis = "PC2",
                      rotatePC = NULL,
                      scale,
                      plotGeno = TRUE,
                      colGeno = NULL,
                      sizeGeno,
                      plotConvHull,
                      plotEnv = TRUE,
                      colorEnvBy,
                      colEnv = NULL,
                      sizeEnv,
                      colorGenoBy,
                      title = NULL) {
  percPC1 <- round(importance[2, primAxis] * 100, 1)
  percPC2 <- round(importance[2, secAxis] * 100, 1)
  if (scale == 1) {
    info <- "environment scaling"
  } else if (scale == 0) {
    info <- "genotype scaling"
  } else if (scale == 0.5) {
    info <- "symmetric scaling"
  } else {
    info <- paste0(round(importance[3, secAxis] * 100, 1), "%")
  }
  if (is.null(title)) {
    title <- paste0(ifelse(GGE, "GGE", "AMMI2"), " biplot for ", trait,
                    " (", info, ") ", year)
  }
  ## Calculate lambda scale.
  lam <- as.numeric(importance[1, c(primAxis, secAxis)])
  lam <- lam * sqrt(nrow(scores))
  lam <- lam ^ scale
  if (plotGeno) {
    ## Create data.frames for genotypes.
    genoDat <- as.data.frame(t(t(scores[, c(primAxis, secAxis)]) / lam))
    genoDat[["type"]] <- "geno"
    if (!is.null(colorGenoBy)) {
      ## Merge the colorGenoBy column to the data.
      genoDat <- merge(genoDat, unique(dat[c("genotype", colorGenoBy)]),
                       by.x = "row.names", by.y = "genotype")
      if (!is.factor(genoDat[[colorGenoBy]])) {
        genoDat[[colorGenoBy]] <- as.factor(genoDat[[colorGenoBy]])
      }
      genoDat <- genoDat[order(genoDat[[colorGenoBy]]), ]
      nColGeno <- nlevels(genoDat[[colorGenoBy]])
      if (length(colGeno) == 0) {
        ## Defaults to black for one color for genotypes.
        ## For more than one colors from statgen.genoColors are used.
        ## Fall back to topo.colors if number of colors in option is too small.
        if (nColGeno == 1) {
          colGeno <- "black"
        } else if (length(getOption("statgen.genoColors")) >= nColGeno) {
          colGeno <- getOption("statgen.genoColors")[1:nColGeno]
        } else {
          colGeno <- topo.colors(nColGeno)
        }
      } else {
        nColGenoArg <- length(colGeno)
        if (nColGenoArg != nColGeno) {
          stop("Number of colors provided doesn't match number of genotype groups:",
               "\n", nColGenoArg, " colors provided, ", nColGeno,
               " groups in data.\n")
        }
      }
      ## Add group variable with contents of colorGenoBy.
      genoDat[[".group"]] <- genoDat[[colorGenoBy]]
      ## Set colors as levels of colorGenoBy.
      levels(genoDat[[colorGenoBy]]) <- colGeno
      rownames(genoDat) <- genoDat[["Row.names"]]
      colnames(genoDat)[colnames(genoDat) == colorGenoBy] <- ".color"
    } else {
      genoDat[[".group"]] <- "genoGroup1"
      genoDat[[".color"]] <- factor(if (is.null(colGeno)) "black" else colGeno[1])
    }
    ## Restrict columns so rbinding to envDat is possible.
    genoDat <- genoDat[c("type", primAxis, secAxis, ".group", ".color")]
    ## Add size, vjust and hjust.
    genoDat[[".size"]] <- sizeGeno
    genoDat[[".vjust"]] <- 1
    genoDat[[".hjust"]] <- 0
  } else {
    genoDat <- NULL
  }
  if (plotEnv) {
    ## Create data.frames for environments.
    envDat <- as.data.frame(t(t(loadings[, c(primAxis, secAxis)]) * lam))
    envDat[["type"]] <- "env"
    if (!is.null(colorEnvBy)) {
      envDat <- merge(envDat, unique(dat[c("trial", colorEnvBy)]),
                      by.x = "row.names", by.y = "trial")
      if (!is.factor(envDat[[colorEnvBy]])) {
        envDat[[colorEnvBy]] <- as.factor(envDat[[colorEnvBy]])
      }
      envDat <- envDat[order(envDat[[colorEnvBy]]), ]
      nColEnv <- nlevels(envDat[[colorEnvBy]])
      if (length(colEnv) == 0) {
        ## Defaults to red for one color for trials.
        ## For more than one colors from statgen.trialColors are used.
        ## Fall back to topo.colors if number of colors in option is too small.
        if (nColEnv == 1) {
          colEnv <- "red"
        } else if (length(getOption("statgen.trialColors")) >= nColEnv) {
          colEnv <- getOption("statgen.trialColors")[1:nColEnv]
        } else {
          colEnv <- topo.colors(nColEnv)
        }
      } else {
        nColEnvArg <- length(colEnv)
        if (nColEnvArg != nColEnv) {
          stop("Number of colors provided doesn't match number of environment ",
               "groups:\n", nColEnvArg, " colors provided, ", nColEnv,
               " groups in data.\n")
        }
      }
      ## Add group variable with contents of colorGenoBy.
      envDat[[".group"]] <- envDat[[colorEnvBy]]
      ## Set colors as levels of colorEnvBy
      levels(envDat[[colorEnvBy]]) <- colEnv
      rownames(envDat) <- envDat[["Row.names"]]
      colnames(envDat)[colnames(envDat) == colorEnvBy] <- ".color"
    } else {
      envDat[[".group"]] <- "envGroup1"
      envDat[[".color"]] <- factor(if (is.null(colEnv)) "red" else colEnv[1])
    }
    ## Restrict columns so rbinding to genoDat is possible.
    envDat <- envDat[c("type", primAxis, secAxis, ".group", ".color")]
    ## Add size, vjust and hjust.
    envDat[[".size"]] <- sizeEnv
    envDat[[".vjust"]] <- "outward"
    envDat[[".hjust"]] <- "outward"
  } else {
    envDat <- NULL
  }
  ## Create total data,frame for convenience.
  totDat <- rbind(genoDat, envDat)
  ## Rotate in such a way that the selected environment or genotype
  ## is aligned with the positive x-axis.
  if (!is.null(rotatePC)) {
    ## Get the angle between the selected env/geno and the positive x-axis.
    xRot <- totDat[rotatePC, primAxis]
    yRot <- totDat[rotatePC, secAxis]
    theta <- atan2(yRot, xRot)
  } else if (GGE) {
    ## Rotation to genotypic main effect as default for GGE.
    xRot <- mean(envDat[, primAxis])
    yRot <- mean(envDat[, secAxis])
    theta <- atan2(yRot, xRot)
  } else {
    ## No rotation by default. Set angle to 0.
    theta <- 0
  }
  ## Rotate clockwise over this angle.
  totDat <- rotatePC(dat = totDat, theta = theta, primAxis = primAxis,
                     secAxis = secAxis)
  genoDat <- rotatePC(dat = genoDat, theta = theta, primAxis = primAxis,
                      secAxis = secAxis)
  envDat <- rotatePC(dat = envDat, theta = theta, primAxis = primAxis,
                     secAxis = secAxis)
  ## Bind together so everything can be plotted in one go.
  ## This has to be done because only one color legend is allowed.
  ## Split between points and text data.
  if (sizeGeno == 0) {
    pointDat <- genoDat
    textDat <- envDat
  } else {
    pointDat <- NULL
    textDat <- rbind(genoDat, envDat)
  }
  ## Get x and y limits and compute plotRatio from them.
  ## This assures 1 unit on the x-asis is identical to 1 unit on the y-axis.
  yMin <- min(totDat[[secAxis]])
  yMax <- max(totDat[[secAxis]])
  xMin <- min(totDat[[primAxis]])
  xMax <- max(totDat[[primAxis]])
  plotRatio <- (xMax - xMin) / (yMax - yMin)
  colGroups <- setdiff(unique(totDat[[".group"]]),
                       c("envGroup1", "genoGroup1"))
  colNamed <- setNames(c(unique(as.character(genoDat[[".color"]])),
                         unique(as.character(envDat[[".color"]]))),
                       unique(as.character(totDat[[".group"]])))
  if (plotGeno && !is.null(colorGenoBy) && plotEnv && !is.null(colorEnvBy)) {
    nGenoGroups <- length(unique(genoDat[[".group"]]))
    nEnvGroups <- length(unique(envDat[[".group"]]))
    nrowGuide <- nGenoGroups
    shapesGuide <- c(rep(16, times = nGenoGroups), rep(NA, times = nEnvGroups))
    sizesGuide <- c(rep(if (sizeGeno == 0) 2 else sizeGeno,
                        times = nGenoGroups),
                    rep(sizeEnv, times = nEnvGroups))
  } else {
    nrowGuide <- shapesGuide <- sizesGuide <- NULL
  }
  lims <- c(min(xMin, yMin), max(xMax, yMax))
  p <- ggplot2::ggplot(totDat, ggplot2::aes(x = .data[[primAxis]],
                                            y = .data[[secAxis]],
                                            color = .data[[".group"]])) +
    ## Needed for a square plot output.
    ggplot2::coord_equal(xlim = lims, ylim = lims, clip = "off") +
    ggplot2::scale_color_manual(breaks = colGroups, values = colNamed,
                                name = NULL,
                                guide = ggplot2::guide_legend(nrow = nrowGuide,
                                                              override.aes = list(shape = shapesGuide,
                                                                                  size = sizesGuide))) +
    ## Add labeling.
    ggplot2::labs(x = paste0(primAxis, " (", percPC1, "%)"),
                  y = paste0(secAxis, " (", percPC2, "%)"),
                  title = title) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   panel.grid = ggplot2::element_blank())
  if (!is.null(pointDat)) {
    ## Plot genotypes as points.
    p <- p + ggplot2::geom_point(data = pointDat,
                                 show.legend = !is.null(colorGenoBy))
  }
  if (!is.null(textDat)) {
    ## Plot genotypes and environments as text.
    p <- p + ggplot2::geom_text(data = textDat,
                                ggplot2::aes(label = rownames(textDat),
                                             size = .data[[".size"]],
                                             vjust = .data[[".vjust"]],
                                             hjust = .data[[".hjust"]]),
                                show.legend = !is.null(colorEnvBy)) +
      ggplot2::scale_size(range = range(textDat[[".size"]]), guide = "none")
  }
  if (plotConvHull) {
    ## Compute convex hull for the points.
    convHulls <- genoDat[grDevices::chull(genoDat[, c(primAxis, secAxis)]), ]
    ## Add convex hull to plot as a polygon.
    p <- p + ggplot2::geom_polygon(color = "darkolivegreen3",
                                   data = convHulls, alpha = 0.2)
    ## Lines from origin perpendicular to convex hull only for GGE2 plots.
    if (GGE) {
      ## Extract x and y coordinates for points on hull. Add first item to the
      ## end to include all edges.
      xConv <- convHulls[[primAxis]]
      yConv <- convHulls[[secAxis]]
      xConv <- c(xConv, xConv[1])
      yConv <- c(yConv, yConv[1])
      ## Compute slopes per segment of the hull.
      slopesConv <- diff(yConv) / diff(xConv)
      ## Compute slopes for the lines perpendicular to the hull.
      slopesPerp <- -1 / slopesConv
      ## Compute the coordinates of the points on the hull through which the
      ## perpendicular lines should go.
      origConv <- yConv[-1] - slopesConv * xConv[-1]
      xNew <- -origConv / (slopesConv - slopesPerp)
      yNew <- slopesPerp * xNew
      ## Expand the lines outward of the hull.
      ## Expansion is done in two steps. First in the x-direction with
      ## computation of the y-coordinate. If this coordinate is outside the
      ## plot area expansion is repeated but then in the y-direction.
      for (i in seq_along(xNew)) {
        if (xNew[i] > 0) {
          xNewI <- lims[2]
          yNewI <- slopesPerp[i] * lims[2]
        } else {
          xNewI <- lims[1]
          yNewI <- slopesPerp[i] * lims[1]
        }
        if (yNewI < lims[1]) {
          yNewI <- lims[1]
          xNewI <- lims[1] / slopesPerp[i]
        } else if (yNewI > lims[2]) {
          yNewI <- lims[2]
          xNewI <- lims[2] / slopesPerp[i]
        }
        xNew[i] <- xNewI
        yNew[i] <- yNewI
      }
      ## Put data for perpendicular lines in a single data set
      ## for ease of plotting.
      perpDat <- data.frame(xend = xNew, yend = yNew)
      ## Add perpendicular lines to plot as segments.
      p <- p + ggplot2::geom_segment(ggplot2::aes(x = 0, y = 0,
                                                  xend = .data[["xend"]],
                                                  yend = .data[["yend"]]),
                                     data = perpDat, col = "grey50",
                                     linewidth = 0.6)
    }
  }
  if (plotEnv) {
    ## Add arrows from origin to environments.
    ## Adding alpha = for transparency causes the arrows not being plotted
    ## after turning off clipping which is needed since labels may fall off
    ## the plot otherwise.
    p <- p +
      ggplot2::geom_segment(data = envDat,
                            ggplot2::aes(x = 0, y = 0,
                                         xend = .data[[primAxis]],
                                         yend = .data[[secAxis]]),
                            arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm")),
                            show.legend = FALSE)
  }
  return(p)
}

#' Helper function for clockwise rotation over an angle of theta.
#'
#' @noRd
#' @keywords internal
rotatePC <- function(dat,
                     primAxis,
                     secAxis,
                     theta) {
  if (!is.null(dat)) {
    rot1 <- dat[[primAxis]] * cos(theta) + dat[[secAxis]] * sin(theta)
    rot2 <- - dat[[primAxis]] * sin(theta) + dat[[secAxis]] * cos(theta)
    dat[[primAxis]] <- rot1
    dat[[secAxis]] <- rot2
  }
  return(dat)
}

#' Extract fitted values.
#'
#' Extract the fitted values for an object of class AMMI.
#'
#' @param object An object of class AMMI
#' @param ... Not used.
#'
#' @return A data.frame with fitted values.
#'
#' @examples
#' ## Run AMMI analysis on TDMaize.
#' geAmmi <- gxeAmmi(TD = TDMaize, trait = "yld")
#'
#' ## Extract fitted values.
#' fitAmmi <- fitted(geAmmi)
#' head(fitAmmi)
#'
#' @family AMMI
#'
#' @export
fitted.AMMI <- function(object,
                        ...) {
  return(object$fitted)
}

#' Extract residuals.
#'
#' Extract the residuals for the fitted AMMI model.
#'
#' @param object An object of class AMMI
#' @param ... Not used.
#'
#' @return A data.frame with residuals.
#'
#' @examples
#' ## Run AMMI analysis on TDMaize.
#' geAmmi <- gxeAmmi(TD = TDMaize, trait = "yld")
#'
#' ## Extract residuals.
#' residAmmi <- residuals(geAmmi)
#' head(residAmmi)
#'
#' @family AMMI
#'
#' @export
residuals.AMMI <- function(object,
                           ...) {
  trait <- object$trait
  fittedGeno <- object$fitted
  residGeno <- merge(fittedGeno, object$dat, by = c("trial", "genotype"))
  residGeno[["residual"]] <- residGeno[["fittedValue"]] - residGeno[[trait]]
  residGeno <- residGeno[c("trial", "genotype", "residual")]
  return(residGeno)
}

#' Report method for class AMMI
#'
#' A pdf report will be created containing a summary of an AMMI object.
#' Simultaneously the same report will be created as a tex file.
#'
#' @param x An object of class AMMI.
#' @param ... Not used.
#' @param outfile A character string, the name and location of the output .pdf
#' and .tex file for the report. If \code{NULL}, a report with a default name
#' will be created in the current working directory.
#'
#' @return A pdf and tex report.
#'
#' @examples
#' ## Run AMMI analysis on TDMaize.
#' geAmmi <- gxeAmmi(TD = TDMaize, trait = "yld")
#' \donttest{
#' ## Create a pdf report summarizing the results.
#' report(geAmmi, outfile = tempfile(fileext = ".pdf"))
#' }
#'
#' @family AMMI
#'
#' @export
report.AMMI <- function(x,
                        ...,
                        outfile = NULL) {
  ## Checks.
  if (nchar(Sys.which("pdflatex")) == 0) {
    stop("An installation of LaTeX is required to create a pdf report.\n")
  }
  createReport(x = x, reportName = "ammiReport.Rnw",
               reportPackage = "statgenGxE", outfile = outfile, ...)
}



