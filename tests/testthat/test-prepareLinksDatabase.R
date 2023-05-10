annot_path <- system.file("exdata/dmel_reference_annotation.gtf.gz", package="LATER")
#annot_path <- "../data/dm6.annot.gtf.gz"
ref_annot <- rtracklayer::import.gff(annot_path)
# test junction reference generation
#junction_path <- "../data/short_read_junctions.SJ.out.tab"
test_that('IsoformDatabase correctly created', {
  expect_error(LATER::prepareIsoformDatabase( ref_annot, tss.window = 50, tes.window = 150),regexp=NA)
})


