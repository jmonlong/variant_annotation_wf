library(GenomicRanges)

load('annotation_database.RData')

region = GRanges("chr7:97357275-104567465")

dgva = subsetByOverlaps(dgva, region)
clinsv = subsetByOverlaps(clinsv, region)
dis.enh = subsetByOverlaps(dis.enh, region)
dis.enh.enc = subsetByOverlaps(dis.enh.enc, region)
sr = subsetByOverlaps(sr, region)

sv.catalogs2 = lapply(sv.catalogs, function(cat){
  subsetByOverlaps(cat, region)
})
names(sv.catalogs2) = names(sv.catalogs)
sv.catalogs = sv.catalogs2

save(dgva, clinsv, dis.enh, dis.enh.enc, sr, sv.catalogs, file='sv_annotation_database.test_subset.RData')
