#' @export
getNews_prBARgraph <- function(...) {
  newsfile <- file.path(system.file(package = "crrBAR"), "NEWS")
  file.show(newsfile)
}
