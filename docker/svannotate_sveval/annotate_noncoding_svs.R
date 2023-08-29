args = commandArgs(TRUE)
## 1: input VCF
## 2: SV database (.RData file)
## 3: ouptput VCF (unzipped VCF)

library(dplyr)
library(sveval)
library(VariantAnnotation)

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
  af.l = info(res.m)$AF
  alt.l = alt.l[which(alt.l>1)]
  res.m = res.m[rep(1:length(alt.l), alt.l)]
  ll = lapply(unlist(lapply(rr.gr$ALT, as.character)), function(x) DNAStringSet(x))
  fixed(res.m)["ALT"] <- DNAStringSetList(ll)
  ## also splits the AF field appropriately
  if(!is.null(af.l)){
    info(res.m)$AF = as(unlist(af.l), "NumericList")
  }
  ## merge back with records that were biallelic
  rbind(res.b, res.m)
}

## add some SV-specific fields that sveval will use
## assumes biallelic VCF object as input
parseSVs <- function(vcf.o){
  rr.gr = rowRanges(vcf.o)
  size.d = unlist(nchar(rr.gr$ALT)) - nchar(rr.gr$REF)
  rowRanges(vcf.o)$size = abs(size.d)
  rowRanges(vcf.o)$type = ifelse(size.d < 0, 'DEL', 'INS')
  rowRanges(vcf.o)$ac = 1
  return(vcf.o)
}

## annotate SVs
annotateFrequency <- function(svs, sv.catalogs, simprep=NULL,
                              min.ol=.1, max.ins.dist=100){
  ## matrix with frequency of each variant (column) in the different SV catalogs (rows)
  af.mat = lapply(sv.catalogs, function(sv.cat){
    af = rep(0, length(svs))
    ol.af = svOverlap(svs, sv.cat, min.ol=min.ol, simprep=simprep,
                      max.ins.dist=max.ins.dist)
    if(is.null(ol.af)){
      return(af)
    }
    ol.af = ol.af %>% 
      as.data.frame %>%
      mutate(af=sv.cat$af[subjectHits]) %>%
      group_by(queryHits) %>% summarize(af=min(af))
    af[ol.af$queryHits] = ol.af$af
    return(af)
  })
  af.mat = do.call(cbind, af.mat)
  ## report the maximum AF across the catalogs
  return(apply(af.mat, 1, max))
}

annotateOverlap <- function(svs, svdb){
  svdb.ol = rep(NA, length(svs))
  ol = svOverlap(svs, svdb)
  if(length(ol)>0){
    ol = ol %>% as.data.frame %>%
      group_by(queryHits) %>% arrange(desc(olScore)) %>% do(head(.,1)) %>%
      dplyr::select(queryHits, subjectHits, olScore)
    svdb.ol[ol$queryHits] = ol$olScore
  }
  svdb.ol
}

annotateSVs <-function(vcf.o){
  ## make a GRanges object
  svs.gr = rowRanges(vcf.o)
  svs.gr$ii = 1:length(vcf.o)
  names(svs.gr) = NULL
  ## ## remove ALT contigs
  chrs=c(1:22, 'X', 'Y')
  ## guess if chromosome names use the 'chr' prefix
  if(all(!(as.character(seqlevels(svs.gr)) %in% chrs)) &
     any(as.character(seqlevels(svs.gr)) %in% paste0('chr', chrs))){
    chrs = paste0('chr', chrs)
  }
  ## keep only non-ALT contigs
  svs.gr = svs.gr[which(as.character(seqnames(svs.gr)) %in% chrs)]
  ## rename tandem duplications into DUP
  svs.gr$type = ifelse(svs.gr$type=='DUP:TANDEM', "DUP", svs.gr$type)
  svs.gr$type = ifelse(svs.gr$type=='DUP', "INS", svs.gr$type)
  ## keep only svs near enhancer of disease genes
  svs.gr = subsetByOverlaps(svs.gr, c(dis.enh, dis.enh.enc), maxgap=10)
  ## annotate frequency and presence in database
  svs.gr$af = suppressWarnings(annotateFrequency(svs.gr, sv.catalogs, simprep=sr))
  svs.gr$clinsv = suppressWarnings(annotateOverlap(svs.gr, clinsv))
  svs.gr$dgv = suppressWarnings(annotateOverlap(svs.gr, dgva))
  ## keep only rare svs
  svs.gr = subset(svs.gr, af<.01 & dgv < .9) 
  ## header for the new INFO fields
  info.h = S4Vectors::DataFrame(
                        Number=rep('1', 3), Type=rep('Float', 3),
                        Description=c(
                          'Allele frequency',
                          'Overlap score with the most similar SV in dbVar Clinical SVs',
                          'Overlap score with the most similar SV in DGV'))
  rownames(info.h) = c('AF', 'CLINSV', 'DGV')
  info(header(vcf.o)) = rbind(info(header(vcf.o)), info.h)
  ## add information back to the vcf object
  vcf.o = vcf.o[svs.gr$ii]
  info(vcf.o)$AF = signif(svs.gr$af, 4)
  info(vcf.o)$CLINSV = NA
  info(vcf.o)$CLINSV = round(svs.gr$clinsv, 4)
  info(vcf.o)$DGV = NA
  info(vcf.o)$DGV = round(svs.gr$dgv, 4)
  ## return new vcf object
  return(vcf.o)
}

##
## 
##

load(args[2])

## read input VCF
vcf.o = readVcf(args[1])

## prepare and annotate SVs
vcf.o = splitMultiAls(vcf.o)
vcf.o = parseSVs(vcf.o)
vcf.o = annotateSVs(vcf.o)

## write annotated VCF
suppressWarnings(writeVcf(vcf.o, file=args[3]))
