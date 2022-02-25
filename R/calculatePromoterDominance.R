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
  annotPairsExp <-
    left_join(
      countData %>% dplyr::select(!gene_id),
      pairsDataBase %>% dplyr::distinct(pairs_id, .keep_all = TRUE),
      by = "pairs_id"
    )
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

