.onLoad <- function(libname = find.package("statgenGxE"),
                    pkgname = "statgenGxE"){
  # CRAN Note avoidance
  if (getRversion() >= "2.15.1")
    utils::globalVariables(c(".", "..count.."))
  invisible()
}
