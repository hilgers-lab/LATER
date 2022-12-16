#' estimateTranscriptionalBias
#' transcriptional bias estimates as in (Alfonso-Gonzalez, C, et al., 2022)
#' @param promoterDominanceTable Table of promoter dominance summaring transcript pairs expression
#' @param method option "chisq" or "fisher" test
#'
#' @return table of significantly biased genes
#' @export
#'
#' @examples
estimateTranscriptionalBias  <-
  function(promoterDominanceTable, method) {
    apa.genes <-
      promoterDominanceTable %>% distinct(tes_id, .keep_all = TRUE) %>% group_by(gene_id) %>%
      filter(n() > 1)  %>% pull(gene_id)
    atss.genes <-
      promoterDominanceTable %>% distinct(promoter_id, .keep_all = TRUE) %>% group_by(gene_id) %>%
      filter(n() > 1)  %>% pull(gene_id)
    atss.apa.genes <- intersect(apa.genes, atss.genes)
    headData <-
      promoterDominanceTable %>% filter(gene_id %in% atss.apa.genes) %>% dplyr::select(gene_id, tes_id, pairs_sum, promoter_id)
    perGeneList <- split(headData, f = headData$gene_id)
    couplingsmatrix <- lapply(perGeneList, function(x) {
     # Produce combination matrix of every TES and Promoter of the gene 
      x1 <-
        x %>% maditr::dcast(tes_id ~ promoter_id, value.var = "pairs_sum") %>% as.data.frame(.)
      rownames(x1) <- x1$tes_id
      x1$tes_id <- NULL
      return(x1)
    })
    couplingsmatrix <-
      lapply(couplingsmatrix, function(x) {
        x[is.na(x)] <- 0
        x
      })
    # Add pseudocount  
    couplingsmatrix <- lapply(couplingsmatrix, function(x) {
      x + 0.7
    })
    # For every per gene matrix perform chisq test and monte carlo simulation to obtain better p-value estimates. 
    if(method=="chisq"){
      perGeneChisqTest <- lapply(couplingsmatrix, function(x) {
        x2 <- chisq.test(x, simulate.p.value = TRUE)
        x1 <- x2$p.value
        return(x1)
      })
      affectedGenes <-
        as.data.frame(do.call(rbind, perGeneChisqTest)) %>% dplyr::rename(p.value.chisq =
                                                                            1) %>% dplyr::mutate(gene_id = rownames(.))
      affectedGenes$p.adj.chisq <-
        p.adjust(affectedGenes$p.value.chisq , method = "BH")
      return(affectedGenes)

    }else{
      # For every per gene matrix perform fisher test and monte carlo simulation to obtain better p-value estimates.  
      perGeneFisher <- lapply(couplingsmatrix, function(x) {
        res.fish <- fisher.test(x, simulate.p.value = TRUE) 
        res.fish <- res.fish$p.value
        return(res.fish)
      })
      resFisher <-
        as.data.frame(do.call(rbind, perGeneFisher)) %>% dplyr::rename(p.value.fisher =
                                                                         1) %>% dplyr::mutate(gene_id = rownames(.))
      resFisher$p.adj.fisher <-
        p.adjust(resFisher$p.value.fisher , method = "BH")
      return(resFisher)
    }
  }
