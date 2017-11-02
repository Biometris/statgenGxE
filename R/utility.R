is_valid_variable_name <- function(x, allow_reserved = TRUE, unique = FALSE) 
{
  # Author: Richie Cotton
  # http://4dpiecharts.com/tag/regex/
  
  ok <- rep.int(TRUE, length(x))
 
  #is name too long?
  #max_name_length <- if(getRversion() < "2.13.0") 256L else 10000L
  #ok[nchar(x) > max_name_length] <- FALSE
 
  #is it a reserved variable, i.e.
  #an ellipsis or two dots then a number?
  if(!allow_reserved)
  {
    ok[x == "..."] <- FALSE    
    ok[grepl("^\\.{2}[[:digit:]]+$", x)] <- FALSE  
  }
 
  #are names valid (and maybe unique)
  ok[x != make.names(x, unique = unique)] <- FALSE
 
  ok
}

char2numeric <- function(x, dec="."){
  # This function is to convert a vector of characters to numeric with specified dec
  if (dec == "."){ 
    y <- as.numeric(x)
  }else{
    y <- sapply(strsplit(x, dec, fixed = TRUE), function(z) as.numeric(paste(z, collapse=".")))
  }
  y
}

file_ext <- function (x) 
{
    # From tools package
    pos <- regexpr("\\.([[:alnum:]]+)$", x)
    ifelse(pos > -1L, substring(x, pos + 1L), "")
}

se.var <- function (x, na.rm = FALSE) 
{
    if (is.matrix(x)) 
        apply(x, 2, se.var, na.rm = na.rm)
    else if (is.vector(x)) {
        if (na.rm) 
            x <- x[!is.na(x)]
        n <- length(x)
        m1 <- sum(x)/length(x)
        m2 <- sum(x^2)/length(x)
        m3 <- sum(x^3)/length(x)
        m4 <- sum(x^4)/length(x)
        sqrt((n*(m4 - 4*m1*m3 + 6*m1*m1*m2 - 3*m1^4)/(n-1) - (n*(m2-m1*m1)/(n-1))^2)/n)                
    }
    else if (is.data.frame(x)) 
        sapply(x, se.var, na.rm = na.rm)
    else se.var(as.vector(x), na.rm = na.rm)
}

skewness <- function (x, na.rm = FALSE) 
{
    if (is.matrix(x)) 
        apply(x, 2, skewness, na.rm = na.rm)
    else if (is.vector(x)) {
        if (na.rm) 
            x <- x[!is.na(x)]
        n <- length(x)
        m1 <- sum(x)/length(x)
        m2 <- sum(x^2)/length(x)
        m3 <- sum(x^3)/length(x)
        (m3 - 3*m1*m2 + 2*m1^3)/(m2 - m1*m1)^(3/2)        
    }
    else if (is.data.frame(x)) 
        sapply(x, skewness, na.rm = na.rm)
    else skewness(as.vector(x), na.rm = na.rm)
}

se.skewness <- function (n)
{
    sqrt((6*n*(n-1))/((n-1)*(n+1)*(n+3)))
}

kurtosis <- function (x, na.rm = FALSE) 
{
    if (is.matrix(x)) 
        apply(x, 2, se.var, na.rm = na.rm)
    else if (is.vector(x)) {
        if (na.rm) 
            x <- x[!is.na(x)]
        n <- length(x)
        m1 <- sum(x)/length(x)
        m2 <- sum(x^2)/length(x)
        m3 <- sum(x^3)/length(x)
        m4 <- sum(x^4)/length(x)
        (m4 - 4*m1*m3 + 6*m1*m1*m2 - 3*m1^4)/(m2 - m1*m1)^2 -3           
    }
    else if (is.data.frame(x)) 
        sapply(x, se.var, na.rm = na.rm)
    else se.var(as.vector(x), na.rm = na.rm)
}

se.kurtosis<- function (n)
{
    sqrt((24*n*(n-1)^2)/((n-2)*(n-3)*(n+5)*(n+3)))
}

wald.test <- function(Sigma, b, positions, df = NULL){
## this function is a cut-down version of wald.test in aod package

  w <- length(positions)
## null hypothesis
  H0 <- rep(0, w)
## L matrix
  L <- matrix(rep(0, length(b) * w), ncol = length(b))
  for(i in 1:w)
    L[i, positions[i]] <- 1
  dimnames(L) <- list(paste("L", as.character(seq(NROW(L))), sep = ""), names(b))
## computations
  f <- L %*% b
  V <- Sigma
  mat <- qr.solve(L %*% V %*% t(L))
  stat <- t(f - H0) %*% mat %*% (f - H0)
  p <- 1 - pchisq(stat, df = w)
  if(is.null(df)){
    res <- list(chi2 = c(chi2 = stat, df = w, P = p))
  }else{
    fstat <- stat / nrow(L)
    df1 <- nrow(L); df2 <- df
    res <- list(chi2 = c(chi2 = stat, df = w, P = p),
                Ftest = c(Fstat = fstat, df1 = df1, df2 = df2, P = 1 - pf(fstat, df1, df2)))
  }
  list(Sigma = Sigma, b = b, positions = positions, L = L, result = res, df = df)
}


