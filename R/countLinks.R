#' count full length reads to 5'-'3 isoforms
#'
#' @param alignments genomic ranges object of alignments (samtools bam to bed)
#' @param linksDatabase 5'-3' isoform link database
#'
#' @return
#' @export
#' @import GenomicAlignments GenomicRanges dplyr edgeR S4Vectors
#' @examples
countLinks <- function(alignments, linksDatabase) {
  # make single nt starts
  startsAlignemnts <-  prepareForCountStarts(alignments, 1)
  startAlignments <-
    readTSSassignment(startsAlignemnts, linksDatabase$TSSCoordinate.base)
  startAlignments <-
    startAlignments %>% as.data.frame(.) %>% dplyr::select(name, promoter_id)
  # make single nt ends
  endsAlignemnts <-  prepareForCountEnds(alignments, 1)
  endsAlignemnts <-
    readTESassignment(endsAlignemnts, linksDatabase$TESCoordinate.base)
  endsAlignemnts <-
    endsAlignemnts %>% as.data.frame(.) %>% dplyr::select(name, tes_id)
  # make pairs
  pairsTested <-
    left_join(endsAlignemnts, startAlignments, by = "name") %>% dplyr::filter(gsub("\\:.*", "", tes_id) == gsub("\\:.*", "", promoter_id))
  pairsTested <-
    left_join(
      pairsTested,
      linksDatabase$pairDataBase %>% dplyr::distinct(tes_id, .keep_all = TRUE) %>% dplyr::select(gene_id, tes_id),
      by = "tes_id"
    )  %>%
    mutate(pairs_id = paste0(
      gene_id,
      ":",
      gsub(".*:", "", pairsTested$promoter_id),
      ":",
      gsub(".*:", "", pairsTested$tes_id)
    ))
  # summarize only reads that expand both pairs
  countedPairs <- pairsTested %>% group_by(pairs_id) %>% tally()
  countedPairsFinal <-
    left_join(
      countedPairs,
      pairsTested %>% dplyr::distinct(pairs_id, .keep_all = TRUE) %>% dplyr::select(gene_id, pairs_id),
      by = "pairs_id"
    )
  # normalize expression in CPMs
  dd <- edgeR::DGEList(counts = countedPairsFinal$n)
  dge <- edgeR::calcNormFactors(dd)
  countedPairsFinal$pairs_cpm <- as.numeric (cpm(dge))
  return(countedPairsFinal)
}
#' prepare counts trimming the reads to their most 5' into a 1 nt window for counting
#'
#' @param x alignments
#' @param window window for trimming
#'
#' @return
#' @export
#'
#' @examples
prepareForCountStarts <- function(x, window) {
  alignments <- x
  pos <- alignments[alignments@strand == "+",]
  neg <- alignments[alignments@strand == "-",]
  end(pos) <- start(pos) + window
  start(neg) <- end(neg) - window
  shortstarts <- c(pos, neg)
  return(shortstarts)
}
#' prepare counts trimming the reads to their most 3' into a 1 nt window for counting
#'
#' @param x alignments
#' @param window window for trimming
#'
#' @return
#' @export
#'
#' @examples
prepareForCountEnds <- function(x, window) {
  alignments <- x
  pos <- alignments[alignments@strand == "+",]
  neg <- alignments[alignments@strand == "-",]
  start(pos) <- end(pos) - window
  end(neg) <- start(neg) + window
  shortstarts <- c(pos, neg)
  return(shortstarts)
}
#' Title
#'
#' @param startsAlignemnts  trimmed reads  output from prepareForCountStarts
#' @param TSSCoordinate.base reference database for TSS
#'
#' @return
#' @export
#'
#' @examples
readTSSassignment <- function(startsAlignemnts, TSSCoordinate.base) {
  tssDb <- TSSCoordinate.base
  ovlps <- findOverlaps(startsAlignemnts , tssDb , maxgap = 50)
  promoterIds <- tssDb[subjectHits(ovlps),]$count
  promoterStarts <- start(tssDb[subjectHits(ovlps),])
  startsAlignemnts2 <- startsAlignemnts[queryHits(ovlps),]
  startsAlignemnts2$promoter_id <- promoterIds
  startsAlignemnts2$promoterStarts <- promoterStarts
  startsAlignemnts2$dist2Assignment <-
    abs(start(startsAlignemnts2) - startsAlignemnts2$promoterStarts)
  startsAlignemnts2 <-
    startsAlignemnts2[order(startsAlignemnts2$dist2Assignment, decreasing = FALSE),]
  startsAlignemnts2 <-
    startsAlignemnts2[!duplicated(startsAlignemnts2$name),]
  return(startsAlignemnts2)
}
#' Title
#'
#' @param endsAlignements trimmed reads  output from prepareForCountEnds
#' @param TESCoordinate.base reference database for TES
#'
#' @return
#' @export
#'
#' @examples
readTESassignment <- function(endsAlignements, TESCoordinate.base) {
  tesDb <- TESCoordinate.base
  ovlps <- findOverlaps(endsAlignements , tesDb , maxgap = 150)
  tesIds <- tesDb[subjectHits(ovlps),]$count
  endStarts <- start(tesDb[subjectHits(ovlps),])
  endsAlignemnts2 <- endsAlignements[queryHits(ovlps),]
  endsAlignemnts2$tes_id <- tesIds
  endsAlignemnts2$endStarts <- endStarts
  endsAlignemnts2$dist2Assignment <-
    abs(start(endsAlignemnts2) - endsAlignemnts2$endStarts)
  endsAlignemnts2 <-
    endsAlignemnts2[order(endsAlignemnts2$dist2Assignment, decreasing = FALSE),]
  endsAlignemnts2 <-
    endsAlignemnts2[!duplicated(endsAlignemnts2$name),]
  return(endsAlignemnts2)
}
