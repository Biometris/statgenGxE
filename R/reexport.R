#' @importFrom statgenSSA createTD
#' @export
statgenSSA::createTD

#' @importFrom statgenSSA fitTD
#' @export
statgenSSA::fitTD

#' @importFrom statgenSSA SSAtoTD
#' @export
statgenSSA::SSAtoTD

#' @importFrom utils getFromNamespace
#' @keywords internal
createSSA <- getFromNamespace(x = "createSSA", ns = "statgenSSA")
