#' @export
#' @useDynLib prBARgraph ccd_ridge
preFiltering <- function(dat) {
  not_all_zeros <- apply(dat, 2, function(x) return(sum(!x == 0) > 0))
  dat = dat[, not_all_zeros]
  return(dat)
}