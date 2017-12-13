# custom tryCatch to return result, errors and warnings
# http://stackoverflow.com/a/24569739/2271856
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

## Extended version of asreml.predict
## Asreml has a bug that may throw a warning message:
## Abnormal termination
## Insufficient workspace - (reset workspace or pworkspace arguments)
## This may be avoided by increasing pworkspace, but this doesn't
## always work.
## If this happens pworkspace is increased in 'small' steps.
predictAsreml <- function(model,
                          classify = "genotype",
                          associate = as.formula("~ NULL"),
                          vcov = TRUE,
                          TD) {
  ## Create tempfile to suppress asreml output messages.
  tmp <- tempfile()
  sink(tmp)
  ## Predict using default settings, i.e. pworkspace = 8e6
  modelP <- tryCatchExt(predict(model, classify = classify,
                                vcov = vcov, associate = associate, data = TD))
  pWorkSpace <- 8e6
  ## While there is a warning, increase pWorkSpace and predict again.
  while (!is.null(modelP$warning) && pWorkSpace < 160e6) {
    pWorkSpace <- pWorkSpace + 8e6
    modelP <- tryCatchExt(predict(model, classify = classify,
                                  vcov = vcov, associate = associate, data = TD,
                                  pworkspace = pWorkSpace))
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

file_ext <- function(x)
{
  # From tools package
  pos <- regexpr(pattern = "\\.([[:alnum:]]+)$", text = x)
  ifelse(pos > -1L, substring(x, pos + 1L), "")
}

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
    se <- sqrt((n * (m4 - 4 * m1 * m3 + 6 * m1 * m1 * m2 -
                       3 * m1 ^ 4) / (n - 1) - (n * (m2 - m1 * m1) / (n - 1)) ^ 2) / n)
  } else if (is.data.frame(x)) {
    se <- sapply(X = x, FUN = seVar, na.rm = na.rm)
  } else {
    se <- seVar(x = as.vector(x), na.rm = na.rm)
  }
  return(se)
}

skewness <- function(x,
                     na.rm = FALSE) {
  if (is.matrix(x)) {
    skw <- apply(X = x, MARGIN = 2, FUN = skewness, na.rm = na.rm)
  } else if (is.vector(x)) {
    if (na.rm) {
      x <- x[!is.na(x)]
    }
    n <- length(x)
    m1 <- sum(x) / length(x)
    m2 <- sum(x ^ 2) / length(x)
    m3 <- sum(x ^ 3) / length(x)
    skw <- (m3 - 3 * m1 * m2 + 2 * m1 ^ 3) / (m2 - m1 * m1) ^ (3 / 2)
  } else if (is.data.frame(x)) {
    skw <- sapply(X = x, FUN = skewness, na.rm = na.rm)
  } else {
    skw <- skewness(x = as.vector(x), na.rm = na.rm)
  }
  return(skw)
}

seSkewness <- function(n) {
  return(sqrt((6 * n * (n - 1)) / ((n - 1) * (n + 1) * (n + 3))))
}

kurtosis <- function(x,
                     na.rm = FALSE) {
  if (is.matrix(x)) {
    kurt <- apply(X = x, MARGIN = 2, FUN = seVar, na.rm = na.rm)
  } else if (is.vector(x)) {
    if (na.rm) {
      x <- x[!is.na(x)]
    }
    n <- length(x)
    m1 <- sum(x) / n
    m2 <- sum(x ^ 2) / n
    m3 <- sum(x ^ 3) / n
    m4 <- sum(x ^ 4) / n
    kurt <- (m4 - 4 * m1 * m3 + 6 * m1 * m1 * m2 - 3 * m1 ^ 4) / (m2 - m1 * m1) ^ 2 - 3
  } else if (is.data.frame(x)) {
    kurt <- sapply(X = x, FUN = seVar, na.rm = na.rm)
  } else {
    kurt <- seVar(x = as.vector(x), na.rm = na.rm)
  }
  return(kurt)
}

seKurtosis <- function(n) {
  return(sqrt((24 * n * (n - 1) ^ 2) / ((n - 2) * (n - 3) * (n + 5) * (n + 3))))
}

waldTest <- function(b,
                     sigma,
                     positions,
                     df = NULL) {
  ## this function is a cut-down version of wald.test in aod package
  w <- length(positions)
  ## L matrix
  L <- matrix(rep(x = 0, times = length(b) * w), ncol = length(b))
  for (i in 1:w) {
    L[i, positions[i]] <- 1
  }
  dimnames(L) <- list(paste0("L", as.character(seq(nrow(L)))), names(b))
  ## computations
  f <- L %*% b
  mat <- qr.solve(Matrix::tcrossprod(L %*% sigma, L))
  stat <- crossprod(f, mat %*% f)
  p <- 1 - pchisq(q = stat, df = w)
  if (is.null(df)) {
    res <- list(chi2 = c(chi2 = stat, df = w, P = p))
  } else {
    fstat <- stat / nrow(L)
    df1 <- nrow(L)
    df2 <- df
    res <- list(chi2 = c(chi2 = stat, df = w, P = p),
                Ftest = c(Fstat = fstat, df1 = df1, df2 = df2,
                          P = 1 - pf(fstat, df1, df2)))
  }
  return(list(sigma = sigma, b = b, positions = positions, L = L, result = res, df = df))
}


