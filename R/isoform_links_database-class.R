# Isoform 5'-3' database class
#' S4 class of promoter database.
#'
#' @slot isoform_database data.frame.
#' @slot TESCoordinate.bins GRanges.
#' @slot TESCoordinate.base GRanges.
#' @slot TSSCoordinate.bins GRanges.
#' @slot TSSCoordinate.base GRanges.
#' @importFrom GenomicRanges GRanges dplyr
#' @return
#' @exportClass IsoformDatabase
#'
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

#' @param x A PromoterAnnotation object
#'
#' @describeIn PromoterAnnotation-class Getter for intronRanges
#' @exportMethod intronRanges

setGeneric("showLinks", function(x) standardGeneric("showLinks"))

#' @describeIn PromoterAnnotation-class Getter for intronRanges
#' @aliases intronRanges,PromoterAnnotation-method

setMethod("showLinks", "IsoformDatabase", function(x) x@isoform_database)

#' @describeIn PromoterAnnotation-class Getter for promoterIdMapping
#' @exportMethod promoterIdMapping

setGeneric("TESCoordinate.bins",
           function(x) standardGeneric("TESCoordinate.bins"))

#' @describeIn PromoterAnnotation-class Getter for promoterIdMapping
#' @aliases promoterIdMapping,PromoterAnnotation-method

setMethod("TESCoordinate.bins", "IsoformDatabase",
          function(x) x@TESCoordinate.bins)

#' @describeIn PromoterAnnotation-class Getter for promoterIdMapping
#' @exportMethod promoterIdMapping

setGeneric("TESCoordinate.base",
           function(x) standardGeneric("TESCoordinate.base"))

#' @describeIn PromoterAnnotation-class Getter for promoterIdMapping
#' @aliases promoterIdMapping,PromoterAnnotation-method

setMethod("TESCoordinate.base", "IsoformDatabase",
          function(x) x@TESCoordinate.base)

#' @describeIn PromoterAnnotation-class Getter for promoterIdMapping
#' @exportMethod promoterIdMapping

setGeneric("TSSCoordinate.base",
           function(x) standardGeneric("TSSCoordinate.base"))

#' @describeIn PromoterAnnotation-class Getter for promoterIdMapping
#' @aliases promoterIdMapping,PromoterAnnotation-method

setMethod("TSSCoordinate.base", "IsoformDatabase",
          function(x) x@TSSCoordinate.base)

#' @describeIn PromoterAnnotation-class Getter for promoterIdMapping
#' @exportMethod promoterIdMapping

setGeneric("TSSCoordinate.bins",
           function(x) standardGeneric("TSSCoordinate.bins"))

#' @describeIn PromoterAnnotation-class Getter for promoterIdMapping
#' @aliases promoterIdMapping,PromoterAnnotation-method

setMethod("TSSCoordinate.bins", "IsoformDatabase",
          function(x) x@TSSCoordinate.bins)


#' @describeIn PromoterAnnotation-class Getter for promoterCoordinates
#' @exportMethod promoterCoordinates

setGeneric("TESCoordinate.base",
           function(x) standardGeneric("TESCoordinate.base"))

#' @describeIn PromoterAnnotation-class Getter for promoterCoordinates
#' @aliases promoterCoordinates,PromoterAnnotation-method

setMethod("TESCoordinate.base", "IsoformDatabase",
          function(x) x@TESCoordinate.base)

###############
### Setters ###

#' @param value intronRanges, promoterIdMapping or promoterCoordinates to
#'   be assigned
#'
#' @describeIn PromoterAnnotation-class Setter for intronRanges
#' @exportMethod 'intronRanges<-'
#' @importFrom methods validObject

setGeneric("isoform_database<-",
           function(x, value) standardGeneric("isoform_database<-"))

#' @describeIn PromoterAnnotation-class Setter for intronRanges
#' @aliases 'intronRanges<-',PromoterAnnotation-method

setMethod("isoform_database<-", "IsoformDatabase", function(x, value) {
  x@isoform_database <- value
  validObject(x)
  x
})

#' @describeIn PromoterAnnotation-class Setter for promoterIdMapping
#' @exportMethod 'promoterIdMapping<-'
#' @importFrom methods validObject

setGeneric("TESCoordinate.base<-",
           function(x, value) standardGeneric("TESCoordinate.base<-"))

#' @describeIn PromoterAnnotation-class Setter for promoterIdMapping
#' @aliases 'promoterIdMapping<-',PromoterAnnotation-method

setMethod("TESCoordinate.base<-", "IsoformDatabase", function(x, value) {
  x@TESCoordinate.base <- value
  validObject(x)
  x
})

setGeneric("TESCoordinate.bins<-",
           function(x, value) standardGeneric("TESCoordinate.bins<-"))

#' @describeIn PromoterAnnotation-class Setter for promoterIdMapping
#' @aliases 'promoterIdMapping<-',PromoterAnnotation-method

setMethod("TESCoordinate.bins<-", "IsoformDatabase", function(x, value) {
  x@TESCoordinate.bins <- value
  validObject(x)
  x
})

setGeneric("TSSCoordinate.bins<-",
           function(x, value) standardGeneric("TSSCoordinate.bins<-"))

#' @describeIn PromoterAnnotation-class Setter for promoterIdMapping
#' @aliases 'promoterIdMapping<-',PromoterAnnotation-method

setMethod("TSSCoordinate.bins<-", "IsoformDatabase", function(x, value) {
  x@TSSCoordinate.bins <- value
  validObject(x)
  x
})

#' @describeIn PromoterAnnotation-class Setter for promoterCoordinates
#' @exportMethod 'promoterCoordinates<-'
#' @importFrom methods validObject

setGeneric("TSSCoordinate.base<-",
           function(x, value) standardGeneric("TSSCoordinate.base<-"))

#' @describeIn PromoterAnnotation-class Setter for promoterCoordinates
#' @aliases 'promoterCoordinates<-',PromoterAnnotation-method

setMethod("TSSCoordinate.base<-", "IsoformDatabase", function(x, value) {
  x@TSSCoordinate.base <- value
  validObject(x)
  x
})

## filtering and handling functions
#' @describeIn TESCoordinate.bins
#' @exportMethod get gene coordiantes

setGeneric("addPromoterDatabase",
           function(x, custom_promoter_annotation,
                    reference_annotation,
                    window) standardGeneric("addPromoterDatabase"))

#' @describeIn PromoterAnnotation-class Getter for promoterIdMapping
#' @aliases promoterIdMapping,PromoterAnnotation-method

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
#' @describeIn TESCoordinate.bins
#' @exportMethod get gene coordiantes

setGeneric("add3pEndDatabase",
           function(x, custom_promoter_annotation,
                    reference_annotation,
                    window) standardGeneric("add3pEndDatabase"))

#' @describeIn PromoterAnnotation-class Getter for promoterIdMapping
#' @aliases promoterIdMapping,PromoterAnnotation-method

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





