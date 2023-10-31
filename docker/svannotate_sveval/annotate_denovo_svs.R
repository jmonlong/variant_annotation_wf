args = commandArgs(TRUE)
## 1: child input VCF
## 2: parent1 input VCF
## 3: parent2 input VCF
## 4: SV database (.RData file)
## 5: ouptput VCF (unzipped VCF)

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
parseSVs <- function(vcf.o){
  rr.gr = rowRanges(vcf.o)
  size.d = unlist(nchar(rr.gr$ALT)) - nchar(rr.gr$REF)
  rowRanges(vcf.o)$size = abs(size.d)
  rowRanges(vcf.o)$type = ifelse(size.d < 0, 'DEL', 'INS')
  rowRanges(vcf.o)$ac = 1
  return(vcf.o)
}

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

annotateSVs <-function(vcf.o, m.sv, f.sv){
  ## make a GRanges object
  svs.gr = rowRanges(vcf.o)
  svs.gr$ii = 1:length(vcf.o)
  names(svs.gr) = NULL
  ## remove ALT contigs
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
  ## annotate de novo candidates
  ## overlap proband SVs with both parents
  ol.m = svOverlap(svs.gr, m.sv, simprep=sr, max.ins.dist=500, min.ol=.05)
  ol.f = svOverlap(svs.gr, f.sv, simprep=sr, max.ins.dist=500, min.ol=.05)
  svs.gr$denovo = TRUE
  svs.gr$denovo[ol.m$queryHits] = FALSE
  svs.gr$denovo[ol.f$queryHits] = FALSE
  svs.gr = subset(svs.gr, denovo)
  ## annotate frequency and presence in database
  svs.gr$af = suppressWarnings(annotateFrequency(svs.gr, sv.catalogs, simprep=sr))
  ## header for the new INFO fields
  info.h = S4Vectors::DataFrame(
                        Number=rep('1', 1), Type=rep('Float', 1),
                        Description=c(
                          'Allele frequency'))
  rownames(info.h) = c('AF')
  info(header(vcf.o)) = rbind(info(header(vcf.o)), info.h)
  ## add information back to the vcf object
  vcf.o = vcf.o[svs.gr$ii]
  info(vcf.o)$AF = signif(svs.gr$af, 4)
  ## return new vcf object
  return(vcf.o)
}

## read input VCF
vcf.o = readVcf(args[1])
m.sv = readSVvcf(args[2])
f.sv = readSVvcf(args[3])

## load annotation database, inc. simple repeats "sr"
load(args[4])

## prepare and annotate SVs
vcf.o = splitMultiAls(vcf.o)
vcf.o = parseSVs(vcf.o)

m.sv = subset(m.sv, type != 'BND')
f.sv = subset(f.sv, type != 'BND')

vcf.o = annotateSVs(vcf.o, m.sv, f.sv)

## write annotated VCF
suppressWarnings(writeVcf(vcf.o, file=args[5]))
