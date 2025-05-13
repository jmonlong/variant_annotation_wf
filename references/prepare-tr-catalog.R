library(sveval)
library(GenomicRanges)
library(dplyr)
library(VariantAnnotation)
library(parallel)
library(tidyr)
library(ggplot2)

overlapSR <- function(tsv.f, add.samp.pref=''){
  sv = read.table(tsv.f, as.is=TRUE)
  colnames(sv) = c('seqnames', 'start', 'ref', 'alt', 'sample')
  sv.dd = sv %>% mutate(end=start+ref, size=alt-ref) %>% 
    makeGRangesFromDataFrame(keep.extra.columns=TRUE)
  ol = findOverlaps(sr, sv.dd, maxgap=1) %>% as.data.frame %>%
    mutate(size=sv.dd$size[subjectHits],
           sample=sv.dd$sample[subjectHits],
           wol=width(pintersect(sr[queryHits], sv.dd[subjectHits])),
           size=ifelse(size<0, -wol, size)) %>% 
    group_by(queryHits, sample) %>% summarize(size=sum(size), .groups='drop')
  sr.sv = sr[ol$queryHits]
  sr.sv$size = ol$size
  sr.sv$sample = paste0(add.samp.pref, ol$sample)
  sr.sv
} 

## simple repeats
sr = read.table('~/Documents/genomedata/simpleRepeat.hg38.txt.gz', as.is=TRUE)
sr = GRanges(sr[,2], IRanges(sr[,3], sr[,4]), motif=sr[,17])
sr = reduce(sr)

message('SV catalog from the HPRC v1.0...')
if(!file.exists('hprc-v1.0-mc.grch38.vcfbub.a100k.wave.svs.sum.tsv.gz')){
  message('download from https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.0-mc.grch38.vcfbub.a100k.wave.vcf.gz')
  message('then: bcftools view -i "STRLEN(REF)>30 | MAX(STRLEN(ALT))>30"')
  message('then: zcat hprc-v1.0-mc.grch38.vcfbub.a100k.wave.svs.vcf.gz | python3 prepare_tr_from_vcf.py | gzip > hprc-v1.0-mc.grch38.vcfbub.a100k.wave.svs.sum.tsv.gz')
}

sr.sv.hprc = overlapSR('hprc-v1.0-mc.grch38.vcfbub.a100k.wave.svs.sum.tsv.gz')

message('SV catalog from the HGSVC2...')
if(!file.exists('hgsvc.freeze3.sv.alt.sum.tsv.gz')){
  message('download from https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/release/v1.0/integrated_callset/freeze3.sv.alt.vcf.gz')
  message('then: zcat hgsvc.freeze3.sv.alt.vcf.gz | python3 prepare_tr_from_vcf.py | gzip > hgsvc.freeze3.sv.alt.sum.tsv.gz')
}

sr.sv.hgsvc = overlapSR('hgsvc.freeze3.sv.alt.sum.tsv.gz', add.samp.pref='HGSVC_')

sr.sv = c(sr.sv.hprc, sr.sv.hgsvc)

sr.sum = sr.sv %>% as.data.frame %>%
  dplyr::select(seqnames, start, end, sample, size) %>%
  pivot_wider(names_from=sample, values_from=size, values_fill=0)

write.table(sr.sum, 'hprc-hgsvc-tr-catalog.tsv', quote=FALSE, row.names=FALSE, sep='\t')
