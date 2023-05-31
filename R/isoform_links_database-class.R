# Isoform 5'-3' database class
#' S4 class of promoter database.
#'
#' @slot isoform_database data.frame.
#' @slot TESCoordinate.bins GRanges.
#' @slot TESCoordinate.bas GRanges.
#' @slot TSSCoordinate.bins GRanges.
#' @slot TSSCoordinate.base GRanges.
#' @importFrom GenomicRanges GRanges
#' @return
#' @exportClass IsoformDatabase
#' @import GenomicFeatures GenomicAlignments S4Vectors dplyr
#' @examples
setClass("IsoformDatabase", slots = c(
  isoform_database = "data.frame",
  TESCoordinate.bins = "GRanges",
  TESCoordinate.base = "GRanges",
  TSSCoordinate.bins = "GRanges",
  TSSCoordinate.base = "GRanges"
),
prototype = list(
  isoform_database = data.frame(),
  TESCoordinate.bins = GenomicRanges::GRanges(),
  TESCoordinate.base = GenomicRanges::GRanges(),
  TSSCoordinate.bins = GenomicRanges::GRanges(),
  TSSCoordinate.base = GenomicRanges::GRanges()
)
)

#' Build isoform database object
#'
#' @param isoform_database
#' @param TESCoordinate.bins
#' @param TESCoordinate.base
#' @param TSSCoordinate.bins
#' @param TSSCoordinate.base
#'
#' @return
#' @export
#'
#' @examples
IsoformDatabase <-
  function(isoform_database = data.frame(),
           TESCoordinate.bins = GenomicRanges::GRanges(),
           TESCoordinate.base = GenomicRanges::GRanges(),
           TSSCoordinate.bins = GenomicRanges::GRanges(),
           TSSCoordinate.base = GenomicRanges::GRanges()) {
    new(
      "IsoformDatabase",
      isoform_database = isoform_database,
      TESCoordinate.bins = TESCoordinate.bins,
      TESCoordinate.base = TESCoordinate.base,
      TSSCoordinate.bins = TSSCoordinate.bins,
      TSSCoordinate.base = TSSCoordinate.base
    )
  }

setValidity("IsoformDatabase", function(object) {
  check <- TRUE
  if (is(object@TESCoordinate.base, 'GRanges') == FALSE) {
    check <- FALSE
  }
  if (is(object@isoform_database, 'data.frame') == FALSE) {
    check <- FALSE
  }
  if (is(object@TSSCoordinate.base, 'GRanges') == FALSE) {
    check <- FALSE
  }
  return(check)
})


#' showLinks
#'
#' @param x
#'
#' @return
#' @export showLinks
#'
#' @examples
setGeneric("showLinks", function(x) standardGeneric("showLinks"))


#' showLinks
#'
#' @param IsoformDatabase
#'
#' @return
#' @export showLinks
#'
#' @examples
setMethod("showLinks", "IsoformDatabase", function(x) x@isoform_database)

#' isoform_database<-
#'
#' @param x
#' @param value
#'
#' @return
#' @export
#'
#' @examples
setGeneric("showLinks<-",
           function(x, value) standardGeneric("showLinks<-"))


#' isoform_database
#'
#' @param IsoformDatabase
#'
#' @return
#' @export
#'
#' @examples
setMethod("showLinks<-", "IsoformDatabase", function(x, value) {
  x@isoform_database <- value
  validObject(x)
  x
})

#' TESCoordinate.bins
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
setGeneric("TESCoordinate.bins",
           function(x) standardGeneric("TESCoordinate.bins"))



#' TESCoordinate.base
#'
#' @param IsoformDatabase
#'
#' @return
#' @export
#'
#' @examples
setMethod("TESCoordinate.bins", "IsoformDatabase",
          function(x) x@TESCoordinate.bins)


#' TESCoordinate.base
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
setGeneric("TESCoordinate.base",
           function(x) standardGeneric("TESCoordinate.base"))



#' TESCoordinate.base M
#'
#' @param IsoformDatabase
#'
#' @return
#' @export
#'
#' @examples
setMethod("TESCoordinate.base", "IsoformDatabase",
          function(x) x@TESCoordinate.base)


#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
setGeneric("TSSCoordinate.base",
           function(x) standardGeneric("TSSCoordinate.base"))


#' TSSCoordinate.base
#'
#' @param IsoformDatabase
#'
#' @return
#' @export
#'
#' @examples
setMethod("TSSCoordinate.base", "IsoformDatabase",
          function(x) x@TSSCoordinate.base)


#' TSSCoordinate.bins
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
setGeneric("TSSCoordinate.bins",
           function(x) standardGeneric("TSSCoordinate.bins"))


#' TSSCoordinate.bins
#'
#' @param IsoformDatabase
#'
#' @return
#' @export
#'
#' @examples
setMethod("TSSCoordinate.bins", "IsoformDatabase",
          function(x) x@TSSCoordinate.bins)



#' TESCoordinate.base
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
setGeneric("TESCoordinate.base",
           function(x) standardGeneric("TESCoordinate.base"))


#' TESCoordinate.base
#'
#' @param IsoformDatabase
#'
#' @return
#' @export
#'
#' @examples
setMethod("TESCoordinate.base", "IsoformDatabase",
          function(x) x@TESCoordinate.base)

###############
### Setters ###


#' isoform_database<-
#'
#' @param x
#' @param value
#'
#' @return
#' @export
#'
#' @examples
setGeneric("isoform_database<-",
           function(x, value) standardGeneric("isoform_database<-"))


#' isoform_database
#'
#' @param IsoformDatabase
#'
#' @return
#' @export
#'
#' @examples
setMethod("isoform_database<-", "IsoformDatabase", function(x, value) {
  x@isoform_database <- value
  validObject(x)
  x
})


#' TESCoordinate.base
#'
#' @param x
#' @param value
#'
#' @return
#' @export
#'
#' @examples
setGeneric("TESCoordinate.base<-",
           function(x, value) standardGeneric("TESCoordinate.base<-"))


#' Title
#'
#' @param IsoformDatabase
#'
#' @return
#' @export
#'
#' @examples
setMethod("TESCoordinate.base<-", "IsoformDatabase", function(x, value) {
  x@TESCoordinate.base <- value
  validObject(x)
  x
})

#' TESCoordinate.bins
#'
#' @param x
#' @param value
#'
#' @return
#' @export
#'
#' @examples
setGeneric("TESCoordinate.bins<-",
           function(x, value) standardGeneric("TESCoordinate.bins<-"))

#' @param IsoformDatabase
#'
#' @describeIn
#' @aliases

setMethod("TESCoordinate.bins<-", "IsoformDatabase", function(x, value) {
  x@TESCoordinate.bins <- value
  validObject(x)
  x
})

#' Title
#'
#' @param x
#' @param value
#'
#' @return
#' @export
#'
#' @examples
setGeneric("TSSCoordinate.bins<-",
           function(x, value) standardGeneric("TSSCoordinate.bins<-"))

#' @describeIn
#' @aliases

setMethod("TSSCoordinate.bins<-", "IsoformDatabase", function(x, value) {
  x@TSSCoordinate.bins <- value
  validObject(x)
  x
})

#' @describeIn

#' @importFrom

setGeneric("TSSCoordinate.base<-",
           function(x, value) standardGeneric("TSSCoordinate.base<-"))

#' @describeIn
#' @aliases

setMethod("TSSCoordinate.base<-", "IsoformDatabase", function(x, value) {
  x@TSSCoordinate.base <- value
  validObject(x)
  x
})

## filtering and handling functions
#' @describeIn TESCoordinate.bins


setGeneric("addPromoterDatabase",
           function(x, custom_promoter_annotation,
                    reference_annotation,
                    window) standardGeneric("addPromoterDatabase"))



#' Title
#'
#' @param IsoformDatabase
#'
#' @return
#' @export
#'
#' @examples
setMethod("addPromoterDatabase", "IsoformDatabase", function(x,
                                                             custom_promoter_annotation,
                                                             reference_annotation,
                                                             window) {
    # get only gene regions
    genes_protein_coding <- reference_annotation[reference_annotation$type == "gene" &
                                              reference_annotation$gene_biotype == "protein_coding"]
    # get only exon regions
    exons_coding <- reference_annotation[reference_annotation$type == "exon" &
                                      reference_annotation$gene_biotype == "protein_coding"]
    # subset promoter database for those not found in reference annotation
    message("identifying TSS not in ref.annotation")
    novel_tss <- subsetByOverlaps(custom_promoter_annotation,
                                  TSSCoordinate.base(x),
                                  maxgap = window,
                                  invert = TRUE)
    message(length(novel_tss)," promoters will be added to the reference annotation" )
    # annotate novel tss to gene
    message("annotate TSSs to gene")
    hits_to_gene <- GenomicRanges::findOverlaps(novel_tss,
                                                genes_protein_coding,
                                                maxgap = window)
    novel_tss <- novel_tss[queryHits( hits_to_gene ), ]
    # add columns to match links database from LATER
    novel_tss$value.group <- seq(1,length(novel_tss))
    novel_tss$gene_id  <- genes_protein_coding[subjectHits( hits_to_gene ), ]$gene_id
    novel_tss <- novel_tss %>%
      as.data.frame(.) %>%
      mutate(count = paste0(gene_id, ":nP",
                            sprintf("%02d", sequence(dplyr::n())))) %>%
      makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
    message("Add novel TSSs to reference database")
    TSSCoordinate.base(x) <- c(TSSCoordinate.base(x),
                                              novel_tss)
    return(x)
  }
  )

## filtering and handling functions
#' @param x
#'
#' @param custom_promoter_annotation
#' @param reference_annotation
#' @param window
#' @export
#' @describeIn TESCoordinate.bins
setGeneric("add3pEndDatabase",
           function(x, custom_promoter_annotation,
                    reference_annotation,
                    window) standardGeneric("add3pEndDatabase"))



#' Title
#'
#' @param IsoformDatabase
#'
#' @return
#' @export
#'
#' @examples
setMethod("add3pEndDatabase", "IsoformDatabase", function(x,
                                                             custom_promoter_annotation,
                                                             reference_annotation,
                                                             window) {
  # get only gene regions
  genes_protein_coding <- reference_annotation[reference_annotation$type == "gene" &
                                                 reference_annotation$gene_biotype == "protein_coding"]
  # get only exon regions
  exons_coding <- reference_annotation[reference_annotation$type == "exon" &
                                         reference_annotation$gene_biotype == "protein_coding"]
  # subset promoter database for those not found in reference annotation
  message("identifying TSS not in ref.annotation")
  novel_tss <- subsetByOverlaps(custom_promoter_annotation,
                                TESCoordinate.base(x),
                                maxgap = window,
                                invert = TRUE)
  message(length(novel_tss)," promoters will be added to the reference annotation" )
  # annotate novel tss to gene
  message("annotate TSSs to gene")
  hits_to_gene <- GenomicRanges::findOverlaps(novel_tss,
                                              genes_protein_coding,
                                              maxgap = window)
  novel_tss <- novel_tss[queryHits( hits_to_gene ), ]
  # add columns to match links database from LATER
  novel_tss$value.group <- seq(1,length(novel_tss))
  novel_tss$gene_id  <- genes_protein_coding[subjectHits( hits_to_gene ), ]$gene_id
  novel_tss <- novel_tss %>%
    as.data.frame(.) %>%
    mutate(count = paste0(gene_id, ":nTE",
                          sprintf("%02d", sequence(dplyr::n())))) %>%
    makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
  message("Add novel TSSs to reference database")
  TSSCoordinate.base(x) <- c(TESCoordinate.base(x),
                             novel_tss)
  return(x)
}
)





#' Title
#'
#' @param x
#' @param GENE_ID
#'
#' @return
#' @export
#'
#' @examples
setGeneric("showGene", function(x,GENE_ID) standardGeneric("showGene"))


#' sshowgene
#'
#' @param IsoformDatabase
#'
#' @return IsoformDatabase
#'
#' @examples
#' @export
setMethod("showGene", "IsoformDatabase", function(x, GENE_ID) {
  showLinks(x) <- showLinks(x) %>% filter(gene_id %in% GENE_ID)
  TESCoordinate.bins(x) <- TESCoordinate.bins(x) %>% filter(gene_id %in% GENE_ID)
  TESCoordinate.base(x) <- TESCoordinate.base(x) %>% filter(gene_id %in% GENE_ID)
  TSSCoordinate.bins(x) <- TSSCoordinate.bins(x) %>% filter(gene_id %in% GENE_ID)
  TSSCoordinate.base(x) <- TSSCoordinate.base(x) %>% filter(gene_id %in% GENE_ID)
  validObject(x)
  x
})


