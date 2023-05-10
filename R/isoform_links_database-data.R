#' Sort by strand
#'
#' @param x GenoicRanges object
#'
#' @return
#' @export
#' @import GenomicFeatures GenomicAlignments S4Vectors dplyr
#' @examples
strandSort <- function(x) {
  c(
    GenomicRanges::sort(x[x@strand == "+"], decreasing = FALSE),
    GenomicRanges::sort(x[x@strand == "-"], decreasing = TRUE)
  )
}
