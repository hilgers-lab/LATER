#plotting functions
#' plotTranscriptionalBias
#'
#' @param LATER
#' @param xlim
#' @param ylim
#' @param residual
#' @param showGenes
#'
#' @return
#' @export
#'
#' @examples
plotTranscriptionalBias <- function(LATER, xlim, ylim, residual,showGenes){
  if(showGenes == TRUE){
    top10Genes <- LATER@result %>%
      filter(p.adj.chisq<0.1) %>%
      arrange(desc(x_squared)) %>%
      head(.,30) %>%
      dplyr::pull(gene_id)
    top10Pair <-  LATER@stats %>%
      dplyr::filter(gene_id %in% top10Genes) %>%
      filter(abs(residuals) > residual) %>%
      arrange(desc(observed), desc(abs(residuals))) %>%
      dplyr::distinct(gene_id, .keep_all = TRUE) %>%
      dplyr::pull(pairs_id)
    LATER@stats %>%
      mutate(outlier_shape = as.factor(ifelse(observed>=xlim | expected>=ylim, 1,2)),
             observed = ifelse(observed>=xlim, xlim, observed),
             expected = ifelse(expected>=ylim, ylim, expected),
             dif = ifelse(abs(residuals) > residual , TRUE, FALSE )) %>%
      ggplot(., aes(x=observed,
                    y=expected,
                    shape = outlier_shape,
                    color = dif,
                    label = ifelse(pairs_id %in% top10Pair, gene_id, ""))) +
      #geom_point(size=3,
      #           alpha=0.7) +
      geom_point(aes(order = dif), size = 3,
                 alpha = 0.7) +
      geom_text(hjust = 0, vjust = 0, size = 3, color="black") +
      theme_classic() +
      xlim(c(1,xlim)) +
      ylim(c(1,ylim)) +
      ylab("log10(Expected Frequency)") +
      xlab("log10(Observed Frequency)") +
      scale_color_manual(values= c("grey", "#E30B5D")) +
      theme(axis.text = element_text(size = 14))
  }else{
    LATER@stats %>%
      mutate(outlier_shape = as.factor(ifelse(log10(observed)>=xlim | log10(expected)>=ylim, TRUE,FALSE)),
             observed = ifelse(log10(observed)>=xlim, xlim, log10(observed)),
             expected = ifelse(log10(expected)>=ylim, ylim, log10(expected)),
             residual = ifelse(abs(residuals) >=residual ,
                               TRUE,
                               FALSE )) %>%
      ggplot(., aes(x=observed,
                    y=expected,
                    shape = outlier_shape,
                    color = residual)) +
      geom_point(size=2,
                 alpha=0.5) +
      theme_classic() +
      xlim(c(0.5,xlim+0.1)) +
      ylim(c(0.5,ylim+0.1)) +
      ylab("log10(Expected Frequency)") +
      xlab("log10(Observed Frequency)") +
      scale_color_manual(values= c("#727B89", "#D70040")) +
      theme(axis.text = element_text(size = 14)) +
      scale_shape(guide = "none")
  }
}

# plot gene

#' plotGeneBias
#'
#' @param LATER
#' @param geneID
#'
#' @return
#' @export
#'
#' @examples
plotGeneBias <- function(LATER, geneID) {
  LATER@dominance %>% filter(gene_id == geneID) %>%
    dplyr::select(c( tes_id ,
                     promoter_id,
                     promoterDominance) ) %>%
    ggplot(., aes(fill=promoter_id,
                  y=promoterDominance,
                  x=tes_id)) +
    geom_bar(position="stack", stat="identity") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45,
                                     vjust = 1,
                                     hjust=1)) +
    ylab("Promoter contribution") +
    scale_fill_aaas()
  }


