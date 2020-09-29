#' @importFrom statgenSTA report
#' @export
statgenSTA::report

#' @importFrom utils getFromNamespace
#' @keywords internal
createSTA <- getFromNamespace(x = "createSTA", ns = "statgenSTA")

#' @keywords internal
createReport <- getFromNamespace(x = "createReport", ns = "statgenSTA")

#' @keywords internal
tryCatchExt <- getFromNamespace(x = "tryCatchExt", ns = "statgenSTA")
