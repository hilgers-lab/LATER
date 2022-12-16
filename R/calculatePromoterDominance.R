#' Caculate promoter dominance
#'
#' @param countData counts data from 5'-3' links
#' @param pairsDataBase database of 5'-3' link isoforms
#'
#' @return
#' @export
#'
#' @examples
calculatePromoterDominance <- function(countData, pairsDataBase) {
  # Assign count data to their respective 3'end and 5'ends.  
  annotPairsExp <-
    left_join(
      countData %>% dplyr::select(!gene_id),
      pairsDataBase %>% dplyr::distinct(pairs_id, .keep_all = TRUE),
      by = "pairs_id"
    )
  # Calculate promoter dominance and other estimates. 
  # promoterDominance: fraction of reads within a TSS-3'end site divided by the total number of the associated 3'end in the gene. 
  # endDominance: Fraction of reads within a given TSS that expressed a given a 3'end. 
  # endFraction: Fraction of 3'end expression over the total number of reads of the gene 
  # startFraction: Fraction of TSS expression compare to total gene expression.
  annotPairsExp <-
  annotPairsExp <-
    annotPairsExp %>% dplyr::group_by(tes_id) %>% dplyr::mutate(end_sum = sum(pairs_cpm)) %>%
    group_by(promoter_id) %>% dplyr::mutate(start_sum = sum(pairs_cpm)) %>% dplyr::group_by(pairs_id) %>% dplyr::mutate(pairs_sum = sum(pairs_cpm)) %>%
    group_by(gene_id) %>% mutate(geneMean = sum(pairs_cpm)) %>% dplyr::mutate(
      promoterDominance = pairs_sum / end_sum,
      endDominance = pairs_sum / start_sum,
      endFraction = end_sum / geneMean,
      startFraction = start_sum / geneMean
    )
  annotPairsExp <- annotPairsExp %>% dplyr::rename(readCounts=n) %>% mutate(pairType=ifelse(is.na(gene_id), "novelPair", "known"))
  return(annotPairsExp)
}

