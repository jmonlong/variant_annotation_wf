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

save(sv.catalogs,
     dgva,
     genc, types.ranked, types.labels,
     clinsv, dis.enh, dis.enh.enc, jarvis, 
     gene.pli, gene.phaplo, gene.ptriplo,
     omim.genes, clingen.genes, clingen.dos.genes, gencc.genes,
     sr, homopol, sd, gaps, map,
     file='sv_annotation_database.test_subset.RData')

df = read.table('hprc-hgsvc-tr-catalog.tsv.gz', as.is=TRUE, sep='\t', header=TRUE)
df = df[overlapsAny(makeGRangesFromDataFrame(df), region),]

write.table(df, file='hprc-hgsvc-tr-catalog.test_subset.tsv', sep='\t', quote=FALSE, row.names=FALSE)
