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

file_ext <- function (x)
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
    m1 <- sum(x) / length(x)
    m2 <- sum(x ^ 2) / length(x)
    m3 <- sum(x ^ 3) / length(x)
    m4 <- sum(x ^ 4) / length(x)
    se <- sqrt((n * (m4 - 4* m1 * m3 + 6 * m1 * m1 * m2 -
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
    skw <- (m3 - 3 * m1 * m2 + 2 * m1 ^ 3) / (m2 - m1 * m1) ^ (3/2)
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

kurtosis <- function (x,
                      na.rm = FALSE) {
  if (is.matrix(x)) {
    kurt <- apply(X = x, MARGIN = 2, FUN = seVar, na.rm = na.rm)
  } else if (is.vector(x)) {
    if (na.rm) {
      x <- x[!is.na(x)]
    }
    n <- length(x)
    m1 <- sum(x) / length(x)
    m2 <- sum(x ^ 2) / length(x)
    m3 <- sum(x ^ 3) / length(x)
    m4 <- sum(x ^ 4) / length(x)
    kurt <- (m4 - 4 * m1 * m3 + 6 * m1 * m1 * m2 - 3 * m1 ^ 4) / (m2 - m1 * m1) ^ 2 - 3
  } else if (is.data.frame(x)) {
    kurt <- sapply(X = x, FUN = seVar, na.rm = na.rm)
  } else {
    kurt <- seVar(x = as.vector(x), na.rm = na.rm)
  }
  return(kurt)
}

seKurtosis<- function(n) {
  return(sqrt((24 * n * (n - 1) ^ 2) / ((n - 2) * (n - 3) * (n + 5) * (n + 3))))
}

waldTest <- function(Sigma,
                     b,
                     positions,
                     df = NULL) {
  ## this function is a cut-down version of wald.test in aod package
  w <- length(positions)
  ## null hypothesis
  H0 <- rep(x = 0, times = w)
  ## L matrix
  L <- matrix(rep(x = 0, times = length(b) * w), ncol = length(b))
  for(i in 1:w) {
    L[i, positions[i]] <- 1
  }
  dimnames(L) <- list(paste0("L", as.character(seq(NROW(L)))), names(b))
  ## computations
  f <- L %*% b
  V <- Sigma
  mat <- qr.solve(L %*% V %*% t(L))
  stat <- t(f - H0) %*% mat %*% (f - H0)
  p <- 1 - pchisq(q = stat, df = w)
  if(is.null(df)) {
    res <- list(chi2 = c(chi2 = stat, df = w, P = p))
  } else {
    fstat <- stat / nrow(L)
    df1 <- nrow(L)
    df2 <- df
    res <- list(chi2 = c(chi2 = stat, df = w, P = p),
                Ftest = c(Fstat = fstat, df1 = df1, df2 = df2,
                          P = 1 - pf(fstat, df1, df2)))
  }
  return(list(Sigma = Sigma, b = b, positions = positions, L = L, result = res, df = df))
}


