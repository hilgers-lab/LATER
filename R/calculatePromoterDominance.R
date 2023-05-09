#' Caculate promoter dominance
#'
#' @param countData counts data from 5'-3' links
#' @param pairsDataBase database of 5'-3' link isoforms
#'
#' @return
#' @export
#'
#' @examples
calculatePromoterDominance <- function(LATER, IsoformDatabase) {
  countData <- isoformCounts(LATER)
  pairsDataBase <- showLinks(IsoformDatabase)
  # Assign count data to their respective 3'end and 5'ends.
  annotPairsExp <-
    left_join(
      countData %>% dplyr::select(!gene_id),
      pairsDataBase %>% dplyr::distinct(pairs_id,
                                        .keep_all = TRUE),
      by = "pairs_id"
    )
  # Calculate promoter dominance and other estimates.
  # promoterDominance: fraction of reads within a TSS-3'end site divided by the total number of the associated 3'end in the gene.
  # endDominance: Fraction of reads within a given TSS that expressed a given a 3'end.
  # endFraction: Fraction of 3'end expression over the total number of reads of the gene
  # startFraction: Fraction of TSS expression compare to total gene expression.
  annotPairsExp <-
    annotPairsExp %>%
    dplyr::group_by(tes_id) %>%
    dplyr::mutate(tes_cpm = sum(cpm)) %>% # note this
    group_by(promoter_id) %>%
    dplyr::mutate(tss_cpm = sum(cpm)) %>%
    dplyr::group_by(pairs_id) %>%
    dplyr::mutate(pairs_cpm = sum(cpm)) %>%
    group_by(gene_id) %>%
    mutate(gene_cpm = sum(cpm)) %>%
    dplyr::mutate(
      promoterDominance = pairs_cpm / tes_cpm,
      endDominance = pairs_cpm / tss_cpm,
      endFraction = tes_cpm / gene_cpm,
      startFraction = tss_cpm / gene_cpm
    )
  annotPairsExp <- annotPairsExp %>%
    dplyr::select(-c(cpm)) %>% dplyr::rename(pairs_read_counts = read_counts) %>%
    mutate(pairType=ifelse(is.na(gene_id), "novelPair", "known"))
  return(annotPairsExp)
}





