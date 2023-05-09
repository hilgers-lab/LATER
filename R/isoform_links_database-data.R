#' Sort by strand
#'
#' @param x GenoicRanges object
#'
#' @return
#' @export
#'
#' @examples
strandSort <- function(x) {
  c(
    GenomicRanges::sort(x[x@strand == "+"], decreasing = FALSE),
    GenomicRanges::sort(x[x@strand == "-"], decreasing = TRUE)
  )
}
