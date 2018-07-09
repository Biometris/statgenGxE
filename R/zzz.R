.onLoad <- function(libname = find.package("RAP"), pkgname = "RAP"){
  # CRAN Note avoidance
  if (getRversion() >= "2.15.1")
    utils::globalVariables(c(".", "..count.."))
  invisible()
}
