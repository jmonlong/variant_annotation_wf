args = commandArgs(TRUE)
## 1: input VCF 1
## 2: input VCF 2
## 3: SV database (.RData file)
## 4: ouptput annotated VCF 1 (unzipped VCF)
## 5: ouptput annotated VCF 2 (unzipped VCF)
## 6: sample name for input VCF 1 (optional)
## 7: sample name for input VCF 2 (optional)

## args = c('../../test/svs.vcf.gz', '../../test/svs.vcf.gz', '../../test/sv_annotation_database.test_subset.RData', 'test1.vcf', 'test2.vcf')

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(sveval))
suppressPackageStartupMessages(library(VariantAnnotation))

##
## Functions
##

## input should be biallelic alleles, but just in case
## split potential multi-allelic records
splitMultiAls <- function(vcf.o){
  rr.gr = rowRanges(vcf.o)
  ## number of ALT alleles
  alt.l = unlist(lapply(nchar(rr.gr$ALT), length))
  res.b = vcf.o[which(alt.l==1)]
  res.m = vcf.o[which(alt.l>1)]
  if(length(res.m)==0){
    ## if no multiallelic records, return input vcf object
    return(vcf.o)
  }
  ## split multi-allelic records
  rr.gr = rowRanges(res.m)
  alt.l = alt.l[which(alt.l>1)]
  res.m = res.m[rep(1:length(alt.l), alt.l)]
  ll = lapply(unlist(lapply(rr.gr$ALT, as.character)), function(x) DNAStringSet(x))
  fixed(res.m)["ALT"] <- DNAStringSetList(ll)
  ## also splits the AF field appropriately
  if(any(colnames(info(res.m)) == 'AF')){
    af.l = info(res.m)$AF
    info(res.m)$AF = as(unlist(af.l), "NumericList")
  }
  ## merge back with records that were biallelic
  rbind(res.b, res.m)
}

## add some SV-specific fields that sveval will use
## assumes biallelic VCF object as input
## also filters to non-ALT contigs
parseSVs <- function(vcf.o){
  rr.gr = rowRanges(vcf.o)
  size.d = unlist(nchar(rr.gr$ALT)) - nchar(rr.gr$REF)
  rowRanges(vcf.o)$size = abs(size.d)
  rowRanges(vcf.o)$type = ifelse(size.d < 0, 'DEL', 'INS')
  rowRanges(vcf.o)$ac = 1
  ## remove ALT contigs
  chrs=c(1:22, 'X', 'Y')
  ## guess if chromosome names use the 'chr' prefix
  if(all(!(as.character(seqlevels(rr.gr)) %in% chrs)) &
     any(as.character(seqlevels(rr.gr)) %in% paste0('chr', chrs))){
    chrs = paste0('chr', chrs)
  }
  ## keep only non-ALT contigs
  nonalt.ii = which(as.character(seqnames(rr.gr)) %in% chrs)
  vcf.o = vcf.o[nonalt.ii]
  ## rename tandem duplications into DUP
  rowRanges(vcf.o)$type = ifelse(rowRanges(vcf.o)$type=='DUP:TANDEM', "DUP", rowRanges(vcf.o)$type)
  rowRanges(vcf.o)$type = ifelse(rowRanges(vcf.o)$type=='DUP', "INS", rowRanges(vcf.o)$type)
  ## add a (shortened) ID field
  new.ids = paste0(as.character(seqnames(vcf.o)), '_',
                   start(vcf.o), '_',
                   end(vcf.o), '_',
                   rowRanges(vcf.o)$type, '_',
                   rowRanges(vcf.o)$size)
  if(is.null(names(vcf.o))){
    rowRanges(vcf.o)$ID = new.ids
  } else {
    rowRanges(vcf.o)$ID = names(vcf.o)
  }
  long.ids = which(nchar(rowRanges(vcf.o)$ID)>50)
  if(length(long.ids)>0){
    rowRanges(vcf.o)$ID[long.ids] = new.ids[long.ids]
  }
  return(vcf.o)
}

annotateSVs <-function(vcf.1, vcf.2){
  ## make GRanges objects
  svs.1 = rowRanges(vcf.1)
  svs.1$ii = 1:length(vcf.1)
  svs.2 = rowRanges(vcf.2)
  svs.2$ii = 1:length(vcf.2)
  names(svs.1) = NULL
  names(svs.2) = NULL
  ## overlap SVs
  ol.o = svOverlap(svs.1, svs.2, simprep=sr, max.ins.dist=500, min.ol=.05)
  ## annotate overlap
  ol.1 = ol.o %>% as.data.frame %>% group_by(queryHits) %>%
    summarize(id=paste(unique(svs.2$ID[subjectHits]), collapse=','))
  svs.1$id2 = NA
  svs.1$id2[ol.1$queryHits] = ol.1$id
  ol.2 = ol.o %>% as.data.frame %>% group_by(subjectHits) %>%
    summarize(id=paste(unique(svs.1$ID[queryHits]), collapse=','))
  svs.2$id2 = NA
  svs.2$id2[ol.2$queryHits] = ol.2$id
  ## header for the new INFO fields
  info.h = S4Vectors::DataFrame(Number='1', Type='String',
                                Description='IDs of the SV(s) in the other VCF. NA if not in the other VCF')
  rownames(info.h) = c('ID2')
  info(header(vcf.1)) = rbind(info(header(vcf.1)), info.h)
  info(header(vcf.2)) = rbind(info(header(vcf.2)), info.h)
  ## add information back to the vcf object
  vcf.1 = vcf.1[svs.1$ii]
  info(vcf.1)$ID2 = svs.1$id2
  vcf.2 = vcf.2[svs.2$ii]
  info(vcf.2)$ID2 = svs.2$id2
  ## return new vcf object
  return(list(vcf.1=vcf.1, vcf.2=vcf.2))
}

## read input VCF 1
if(length(args) == 7){
  param <- ScanVcfParam(samples=args[6])
  vcf.1 = readVcf(args[1], param=param)
} else {
  vcf.1 = readVcf(args[1])
}
## prepare and annotate SVs
vcf.1 = splitMultiAls(vcf.1)
vcf.1 = parseSVs(vcf.1)

## read input VCF 2
if(length(args) == 7){
  param <- ScanVcfParam(samples=args[7])
  vcf.2 = readVcf(args[2], param=param)
} else {
  vcf.2 = readVcf(args[2])
}
## prepare and annotate SVs
vcf.2 = splitMultiAls(vcf.2)
vcf.2 = parseSVs(vcf.2)

## load annotation database, inc. simple repeats "sr"
load(args[3])

## overlap SVs
ann.o = annotateSVs(vcf.1, vcf.2)

## write annotated VCF
suppressWarnings(writeVcf(ann.o$vcf.1, file=args[4]))
suppressWarnings(writeVcf(ann.o$vcf.2, file=args[5]))
