#' estimatePromoterDominance
#'
#' @param LATER
#' @param IsoformDatabase
#' @param method
#'
#' @return
#' @export estimatePromoterDominance
estimatePromoterDominance <- function(LATER, IsoformDatabase, method) {
  dominance <- calculatePromoterDominance(LATER, IsoformDatabase)
  transcriptional_bias <- estimateTranscriptionalBias(dominance, method)
  dominance(LATER) <- dominance
  result(LATER) <- transcriptional_bias$affectedGenes
  stats(LATER) <- transcriptional_bias$stats
  return(LATER)
  }
