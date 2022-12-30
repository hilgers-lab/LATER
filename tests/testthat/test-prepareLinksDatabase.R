annot_path <- system.file("exdata/dm6.annot.gtf.gz", package="SaiLoR")
#annot_path <- "../data/dm6.annot.gtf.gz"
ref_annot <- rtracklayer::import.gff(annot_path)
# test junction reference generation
#junction_path <- "../data/short_read_junctions.SJ.out.tab"
test_that('LinkDatabase correctly created', {
  expect_error(prepareLinksDatabase( ref_annot, tss.window = 50, tes.window = 150),regexp=NA)
})


