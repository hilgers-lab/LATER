
#' prepareIsoformDatabase
#' Create a reference database of all TSS-TES link associations based on annotation
#'
#' @param annotation reference annotation GenomicRanges object
#' @param tss.window window to consider distinct TSS
#' @param tes.window window to consider distinct TESassemblies
#'
#' @return promoterDatabase with classification of links
#' @import GenomicFeatures GenomicAlignments S4Vectors dplyr tibble
#' @export prepareIsoformDatabase
prepareIsoformDatabase <- function(annotation, tss.window, tes.window) {
  Isoform_Database <- IsoformDatabase()
  # Build 5'-3' links data base
  # Exon ids
  txdb <- GenomicFeatures::makeTxDbFromGRanges(annotation)
  ebt <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)
  t2g <- AnnotationDbi::select(txdb,
                               keys = names(ebt),
                               keytype = "TXNAME",
                               columns = "GENEID"
  )
  e2 <- BiocGenerics::unlist(ebt)
  e2$transcript_id <- names(e2)
  e2$gene_id <- t2g$GENEID[match(e2$transcript_id, t2g$TXNAME)]
  e2$exon_id <- e2$exon_name
  e2$exon_name <- NULL
  e2$type <- "exon"
  names(e2) <- NULL
  mcols(e2) <- mcols(e2)[, c(
    "exon_id", "exon_rank",
    "transcript_id", "gene_id", "type"
  )]
  bins <- list()
  # TSS data base
  # take first position per transcript and make it single nt
  tss.bins <-
    strandSort(
      plyranges::mutate(
        plyranges::anchor_5p(
          dplyr::filter(e2, exon_rank == 1)), width = 1))
  # make unique TSS starts merging in a 50nt window.
  message("Preparing TSSBase")
  tss.base <-
    strandSort(
      GenomicRanges::makeGRangesFromDataFrame(
        reshape::melt(GenomicRanges::reduce(
          GenomicRanges::split(tss.bins, ~gene_id),
          min.gapwidth = tss.window
        )),
        keep.extra.columns = TRUE
      )
    )
  bins$tss.bins <- tss.bins
  bins$tss.base <- tss.base
  tss.base <-
    tibble::as_tibble(tss.base) %>%
    dplyr::rename(gene_id = value.group_name) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::mutate(
      count = paste0(gene_id, ":P",
                     sprintf("%02d", sequence(dplyr::n())))) %>%
    GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
  # annotate isoforms with promoter_id
  ii <- GenomicRanges::findOverlaps(tss.bins, tss.base, maxgap = tss.window - 1)
  tss.bins.annot <-
    GenomicRanges::makeGRangesFromDataFrame(rbind(data.frame(
      tss.bins[S4Vectors::queryHits(ii)],
      tss.base[S4Vectors::subjectHits(ii)]
    )), keep.extra.columns = TRUE)
  tss.bins.annot <-
    tss.bins.annot %>%
    dplyr::mutate(
      transcript_id = transcript_id,
      promoter_id = count
    ) %>%
    plyranges::select(gene_id, transcript_id, promoter_id)
  # TES data base
  # last exon per transcript
  le <-
    GenomicRanges::makeGRangesFromDataFrame(
      e2 %>%
        as.data.frame(.) %>%
        group_by(transcript_id) %>%
        dplyr::filter(exon_rank == max(exon_rank)),
      keep.extra.columns = TRUE
    )
  tes.bins <- strandSort(
    plyranges::mutate(plyranges::anchor_3p(le), width = 1))
  # make unique TES starts merging in a 50nt window.
  message("Preparing TES database")
  tes.base <-
    strandSort(
      GenomicRanges::makeGRangesFromDataFrame(
        reshape::melt(GenomicRanges::reduce(
          (
            GenomicRanges::split(tes.bins, ~gene_id)
          ),
          min.gapwidth =
            tes.window
        )),
        keep.extra.columns = TRUE
      )
    )
  tes.base <-
    tibble::as.tibble(tes.base) %>%
    dplyr::rename(gene_id = value.group_name) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::mutate(
      count = paste0(gene_id, ":TE",
                     sprintf("%02d", sequence(dplyr::n())))) %>%
    GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
  # assign tes_ids to isoforms
  message("Assign TES to isoforms")
  ii <- findOverlaps(tes.bins, tes.base, maxgap = tes.window - 1)
  tes.bins.annot <-
    GenomicRanges::makeGRangesFromDataFrame(rbind(data.frame(
      tes.bins[S4Vectors::queryHits(ii)],
      tes.base[S4Vectors::subjectHits(ii)]
    )), keep.extra.columns = TRUE)
  tes.bins.annot <-
    tes.bins.annot %>%
    dplyr::mutate(transcript_id = transcript_id,
                  tes_id = count) %>%
    plyranges::select(gene_id, transcript_id, tes_id)
  bins$tes.bins <- tes.bins
  bins$tes.base <- tes.base
  # create link database
  message("Create TSS/TES links database")
  linksDbs <- dplyr::left_join(
    as.data.frame(tes.bins.annot),
    as.data.frame(tss.bins.annot),
    by = "transcript_id"
  ) %>%
    dplyr::select(gene_id.x, transcript_id, promoter_id, tes_id) %>%
    dplyr::rename(gene_id = gene_id.x) %>%
    dplyr::mutate(pairs_id = paste(promoter_id, gsub(".*:", "", tes_id), sep = ":"))
  # Classify promoter type and utr type
  # sorted with proximal first distal later
  # output: class all 3'end isoforms
  le.sort <- strandSort(le)
  tt <- linksDbs %>%
    dplyr::group_by(gene_id) %>%
    dplyr::mutate(
      utr_type = dplyr::case_when(
        tes_id == min(tes_id) ~ "proximal",
        tes_id == max(tes_id) ~ "distal",
        TRUE ~ "other"
      )
    )
  tt <- tt %>%
    group_by(gene_id) %>%
    dplyr::mutate(
      promoter_type = case_when(
        promoter_id == min(promoter_id) ~ "distal",
        promoter_id == max(promoter_id) ~ "proximal",
        TRUE ~ "intermediate"
      )
    )
  # classify APA-ATSS genes
  # genes with more than 1 promoter different promoter
  message("Classify genes")
  atss.gene <- tt %>%
    dplyr::distinct(promoter_id, .keep_all = TRUE) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::pull(gene_id)
  apa.gene <- tt %>%
    dplyr::distinct(tes_id, .keep_all = TRUE) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::pull(gene_id)
  tt <-
    tt %>% dplyr::mutate(
      tss.status = ifelse(gene_id %in% atss.gene, "ATSS", "single_promoter"),
      apa.status = ifelse(gene_id %in% apa.gene, "APA", "noAPA")
    )
  # # # remove neibouring genes missassignments
  tt$tes_gene <- gsub("\\:.*", "", tt$tes_id)
  tt <- subset(tt, gene_id == tes_gene)
  tt$tes_gene <- NULL
  tt$promoter_gene <- gsub("\\:.*", "", tt$promoter_id)
  tt <- subset(tt, gene_id == promoter_gene)
  tt$promoter_gene <- NULL
  bins$links <- tt
  tt <- tt %>%
    dplyr::group_by(gene_id) %>%
    dplyr::mutate(
      utr_type = dplyr::case_when(
        tes_id == min(tes_id) ~ "proximal",
        tes_id == max(tes_id) ~ "distal",
        TRUE ~ "other"
      )
    )
  tt <- tt %>%
    group_by(gene_id) %>%
    dplyr::mutate(
      promoter_type = case_when(
        promoter_id == min(promoter_id) ~ "distal",
        promoter_id == max(promoter_id) ~ "proximal",
        TRUE ~ "intermediate"
      )
    )
  atss.gene <- tt %>%
    dplyr::distinct(promoter_id, .keep_all = TRUE) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::pull(gene_id)
  apa.gene <- tt %>%
    dplyr::distinct(tes_id, .keep_all = TRUE) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::pull(gene_id)
  tt <-
    tt %>% dplyr::mutate(
      tss.status = ifelse(gene_id %in% atss.gene, "ATSS", "single_promoter"),
      apa.status = ifelse(gene_id %in% apa.gene, "APA", "noAPA")
    )
  result <- list()
  result$pairDataBase <- tt
  result$TESCoordinate.bins <- tes.bins
  result$TESCoordinate.base <- tes.base
  result$TSSCoordinate.bins <- tss.bins
  result$TSSCoordinate.base <- tss.base
  message("Database Sucessfully created")
  ### Build Isoform Database
  complete.isoformDatabase <- IsoformDatabase(
    isoform_database = result$pairDataBase,
    TESCoordinate.bins = result$TESCoordinate.bins,
    TESCoordinate.base = result$TESCoordinate.base,
    TSSCoordinate.bins = result$TSSCoordinate.bins,
    TSSCoordinate.base = result$TSSCoordinate.base
  )
  return(complete.isoformDatabase)
}







