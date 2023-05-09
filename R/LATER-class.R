# LATER class
#' S4 class of promoter database.
#'
#' @slot result data.frame.
#' @slot stats data.frame.
#' @slot readAssignments data.frame.
#' @slot isoformCounts data.frame.
#' @importFrom GenomicRanges GRanges dplyr
#' @return
#' @exportClass LATER
#'
#' @examples
setClass("LATER", slots = c(
  result = "data.frame",
  stats = "data.frame",
  dominance = "data.frame",
  readAssignments = "data.frame",
  isoformCounts = "data.frame"
),
prototype = list(
  result = data.frame(),
  stats = data.frame(),
  dominance = data.frame(),
  readAssignments = data.frame(),
  isoformCounts = data.frame()
)
)

#' Builder LATER
#'
#' @param result
#' @param stats
#' @param readAssignments
#' @param isoformCounts
#'
#' @return
#' @export
#'
#' @examples
LATER <-
  function(result = data.frame(),
           stats = data.frame(),
           dominance = data.frame(),
           readAssignments = data.frame(),
           isoformCounts = data.frame()) {
    new(
      "LATER",
      result = result,
      stats = stats,
      dominance = dominance,
      readAssignments = readAssignments,
      isoformCounts = isoformCounts
    )
  }


setValidity("LATER", function(object) {
  check <- TRUE
  if (is(object@result, 'data.frame') == FALSE) {
    check <- FALSE
  }
  if (is(object@stats, 'data.frame') == FALSE) {
    check <- FALSE
  }
  if (is(object@readAssignments, 'data.frame') == FALSE) {
    check <- FALSE
  }
  return(check)
})

#' @param x A PromoterAnnotation object
#'
#' @describeIn PromoterAnnotation-class Getter for intronRanges
#' @exportMethod intronRanges

setGeneric("result", function(x) standardGeneric("result"))

#' @describeIn PromoterAnnotation-class Getter for intronRanges
#' @aliases intronRanges,PromoterAnnotation-method

setMethod("result", "LATER", function(x) x@result)

#' @describeIn PromoterAnnotation-class Getter for promoterIdMapping
#' @exportMethod promoterIdMapping

setGeneric("stats",
           function(x) standardGeneric("stats"))

#' @describeIn PromoterAnnotation-class Getter for promoterIdMapping
#' @aliases promoterIdMapping,PromoterAnnotation-method

setMethod("stats", "LATER",
          function(x) x@stats)

#' @describeIn PromoterAnnotation-class Getter for promoterIdMapping
#' @exportMethod promoterIdMapping

setGeneric("dominance",
           function(x) standardGeneric("dominance"))

#' @describeIn PromoterAnnotation-class Getter for promoterIdMapping
#' @aliases promoterIdMapping,PromoterAnnotation-method

setMethod("dominance", "LATER",
          function(x) x@dominance)

#' @describeIn PromoterAnnotation-class Getter for promoterIdMapping
#' @exportMethod promoterIdMapping

setGeneric("readAssignments",
           function(x) standardGeneric("readAssignments"))

#' @describeIn PromoterAnnotation-class Getter for promoterIdMapping
#' @aliases promoterIdMapping,PromoterAnnotation-method

setMethod("readAssignments", "LATER",
          function(x) x@readAssignments)

#' @describeIn PromoterAnnotation-class Getter for promoterIdMapping
#' @exportMethod promoterIdMapping

setGeneric("isoformCounts",
           function(x) standardGeneric("isoformCounts"))

#' @describeIn PromoterAnnotation-class Getter for promoterIdMapping
#' @aliases promoterIdMapping,PromoterAnnotation-method

setMethod("isoformCounts", "LATER",
          function(x) x@isoformCounts)



###############
### Setters ###

#' @param value intronRanges, promoterIdMapping or promoterCoordinates to
#'   be assigned
#'
#' @describeIn PromoterAnnotation-class Setter for intronRanges
#' @exportMethod 'intronRanges<-'
#' @importFrom methods validObject

setGeneric("result<-",
           function(x, value) standardGeneric("result<-"))

#' @describeIn PromoterAnnotation-class Setter for intronRanges
#' @aliases 'intronRanges<-',PromoterAnnotation-method

setMethod("result<-", "LATER", function(x, value) {
  x@result <- value
  validObject(x)
  x
})

#' @describeIn PromoterAnnotation-class Setter for promoterIdMapping
#' @exportMethod 'promoterIdMapping<-'
#' @importFrom methods validObject

setGeneric("stats<-",
           function(x, value) standardGeneric("stats<-"))

#' @describeIn PromoterAnnotation-class Setter for promoterIdMapping
#' @aliases 'promoterIdMapping<-',PromoterAnnotation-method

setMethod("stats<-", "LATER", function(x, value) {
  x@stats <- value
  validObject(x)
  x
})


setGeneric("dominance<-",
           function(x, value) standardGeneric("dominance<-"))

#' @describeIn PromoterAnnotation-class Setter for intronRanges
#' @aliases 'intronRanges<-',PromoterAnnotation-method

setMethod("dominance<-", "LATER", function(x, value) {
  x@dominance <- value
  validObject(x)
  x
})

setGeneric("readAssignments<-",
           function(x, value) standardGeneric("readAssignments<-"))

#' @describeIn PromoterAnnotation-class Setter for promoterIdMapping
#' @aliases 'promoterIdMapping<-',PromoterAnnotation-method

setMethod("readAssignments<-", "LATER", function(x, value) {
  x@readAssignments <- value
  validObject(x)
  x
})

setGeneric("isoformCounts<-",
           function(x, value) standardGeneric("isoformCounts<-"))

#' @describeIn PromoterAnnotation-class Setter for promoterIdMapping
#' @aliases 'promoterIdMapping<-',PromoterAnnotation-method

setMethod("isoformCounts<-", "LATER", function(x, value) {
  x@isoformCounts <- value
  validObject(x)
  x
})



##
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

