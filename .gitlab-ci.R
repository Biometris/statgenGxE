installIfNeeded <- function(pkg,
                            repos,
                            quiet) {
  pkgPath <- find.package(pkg, quiet = quiet)
  if (length(pkgPath) == 0) {
    message("NOTE: pkg ", pkg, " missing, installing...")
    install.packages(pkg, repos = repos, quiet = quiet)
  }
}

pkgsUpdate <- function(repos = "https://cran.rstudio.com",
                       quiet = TRUE,
                       instPkgdown = FALSE) {
  installIfNeeded(pkg = "remotes", repos = repos, quiet = quiet)
  if (instPkgdown) {
    installIfNeeded(pkg = "pkgdown", repos = repos, quiet = quiet)
  }
  remotes::install_deps(repos = repos, quiet = quiet, dependencies = TRUE)
  cat("INSTALLED:\n")
  instld <- as.data.frame(installed.packages())
  rownames(instld) <- NULL
  print(instld[, c("Package", "Version")])
  return(invisible(TRUE))
}


