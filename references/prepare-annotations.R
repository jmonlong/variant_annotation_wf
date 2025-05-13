## packages to read SVs in VCF and manipulate genomic ranges
library(sveval)
library(GenomicRanges)
library(dplyr)

message('gnomAD-SV...')
if(!file.exists('nstd166.GRCh38.variant_call.vcf.gz')){
  download.file('https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/nstd166.GRCh38.variant_call.vcf.gz',
                'nstd166.GRCh38.variant_call.vcf.gz')
}
gnomad = readSVvcf('nstd166.GRCh38.variant_call.vcf.gz', other.field='AF')
gnomad = subset(gnomad, AF>0)
gnomad.cnv = subset(gnomad, type=='CNV')
gnomad.cnv$type = 'DEL'
gnomad = c(gnomad, gnomad.cnv)
gnomad$type = ifelse(gnomad$type=='CNV', 'DUP', gnomad$type)
gnomad$type = ifelse(gnomad$type=="DUP", "INS", gnomad$type)
gnomad$af = gnomad$AF
gnomad$AF = NULL
seqlevels(gnomad) = paste0('chr', seqlevels(gnomad))
## table(gnomad$type)

message('Gene annotation...')
types.ranked = c('CDS', 'UTR', 'promoter', 'gene')
types.labels = c('coding', 'UTR', 'promoter', 'intronic')
if(!file.exists('gencode.v41.annotation.rds')){
  if(!file.exists('gencode.v41.annotation.gtf.gz')){
    download.file('https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gtf.gz', 'gencode.v41.annotation.gtf.gz')
  }
  genc = rtracklayer::import('gencode.v41.annotation.gtf.gz')
  genc = subset(genc, type %in% types.ranked)
  prom = promoters(subset(genc, type=='gene'))
  prom$type = 'promoter'
  genc = c(genc, prom)
  mcols(genc) = mcols(genc)[,c('type', 'gene_name', 'gene_type')]
  saveRDS(genc, 'gencode.v41.annotation.rds')
} else{
  genc = readRDS('gencode.v41.annotation.rds')
}

message('pli scores')
if(!file.exists('gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz')){
  download.file('https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz', 'gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz')
}
pli.df = read.table('gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz', as.is=TRUE, header=TRUE, sep='\t')
pli.df = pli.df %>% filter(!duplicated(gene))
gene.pli = pli.df$pLI
names(gene.pli) = pli.df$gene

message('Clinical SVs...')
if(!file.exists('nstd102.GRCh38.variant_call.tsv.gz')){
  download.file('https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/tsv/nstd102.GRCh38.variant_call.tsv.gz', 'nstd102.GRCh38.variant_call.tsv.gz')
}
clinsv = read.table('nstd102.GRCh38.variant_call.tsv.gz',
                    header=TRUE, as.is=TRUE, sep='\t', skip=1, comment='', quote='')
gain.types = c('copy number gain', 'duplication', 'insertion', 'tandem duplication')
loss.types = c('copy number loss', 'deletion')
clinsv = clinsv %>% filter(clinical_significance=='Pathogenic',
                           variant_call_type %in% c(gain.types, loss.types),
                           chr %in% c(1:22, 'X','Y')) %>%
  mutate(chr=paste0('chr', chr),
         type=ifelse(variant_call_type %in% gain.types, 'INS', 'DEL'),
         start=ifelse(is.na(start), ifelse(is.na(inner_start),
                                           outer_start, inner_start), start),
         end=ifelse(is.na(stop), ifelse(is.na(inner_stop),
                                        outer_stop, inner_stop), stop),
         size=end-start) %>%
  select(chr, start, end, size, type) %>% makeGRangesFromDataFrame(keep.extra.columns=TRUE)

message("Dosage sensitivity map...")
# from https://zenodo.org/record/6347673
if(!file.exists('Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz')){
  download.file('https://zenodo.org/record/6347673/files/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz?download=1', 'Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz')
}
gene.dosage = read.table('Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz', as.is=TRUE, header=FALSE)
colnames(gene.dosage) = c('gene', 'pHaplo', 'pTriplo')
gene.dosage = subset(gene.dosage, !duplicated(gene))
gene.phaplo = gene.dosage$pHaplo
gene.ptriplo = gene.dosage$pTriplo
names(gene.ptriplo) = names(gene.phaplo) = gene.dosage$gene

message("SV catalog from Giraffe-SV paper...")
if(!file.exists('vggiraffe-sv-2504kgp-svsites.tsv.gz')){
  download.file('https://cgl.gi.ucsc.edu/data/giraffe/calling/vggiraffe-sv-2504kgp-svsites.tsv.gz', 'vggiraffe-sv-2504kgp-svsites.tsv.gz')
}
sv.giraffe = read.table('vggiraffe-sv-2504kgp-svsites.tsv.gz', as.is=TRUE, header=TRUE)
sv.giraffe = sv.giraffe %>% select(seqnames, start, end, type, size, af) %>%
  makeGRangesFromDataFrame(keep.extra.columns=TRUE)
sv.giraffe$type = ifelse(sv.giraffe$type=="DUP", "INS", sv.giraffe$type)

message('SV catalog from the HPRC v1.0...')
if(!file.exists('hprc-v1.0-mc.grch38.vcfbub.a100k.wave.svs.vcf.gz')){
  message('download from https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.0-mc.grch38.vcfbub.a100k.wave.vcf.gz')
  message('then: ')
  message('bcftools view -i "STRLEN(REF)>30 | MAX(STRLEN(ALT))>30" hprc-v1.0-mc.grch38.vcfbub.a100k.wave.vcf.gz | gzip > hprc-v1.0-mc.grch38.vcfbub.a100k.wave.svs.vcf.gz')
}
sv.hprc = readSVvcf.multisamps('hprc-v1.0-mc.grch38.vcfbub.a100k.wave.svs.vcf.gz')
sv.hprc$ref = sv.hprc$alt = NULL
sv.hprc$ac = ifelse(sv.hprc$ac==0, 1, sv.hprc$ac)
af0 = sort(unique(sv.hprc$af))[2]
sv.hprc$af = ifelse(sv.hprc$af==0, af0, sv.hprc$af)
sv.hprc$type = ifelse(sv.hprc$type=="DUP", "INS", sv.hprc$type)

message('SV catalog from the HGSVC2...')
if(!file.exists('hgsvc.freeze3.sv.alt.vcf.gz')){
  download.file('https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/release/v1.0/integrated_callset/freeze3.sv.alt.vcf.gz', 'hgsvc.freeze3.sv.alt.vcf.gz')
}
sv.hgsvc = readSVvcf.multisamps('hgsvc.freeze3.sv.alt.vcf.gz')
sv.hgsvc$ref = sv.hgsvc$alt = NULL

message('SV catalog from the 1000GP ONT Miller...')
if(!file.exists('20240423_jasmine_intrasample_noBND_custom_suppvec_alphanumeric_header_JASMINE.vcf.gz')){
  download.file('https://s3.amazonaws.com/1000g-ont/Gustafson_etal_2024_preprint_SUPPLEMENTAL/20240423_jasmine_intrasample_noBND_custom_suppvec_alphanumeric_header_JASMINE.vcf.gz', '20240423_jasmine_intrasample_noBND_custom_suppvec_alphanumeric_header_JASMINE.vcf.gz')
}
sv.ontkgp = readSVvcf.multisamps('20240423_jasmine_intrasample_noBND_custom_suppvec_alphanumeric_header_JASMINE.vcf.gz')
sv.ontkgp$af = sv.ontkgp$ac / (2* (sv.ontkgp$nrefs + sv.ontkgp$ncalls))
sv.ontkgp$ref = sv.ontkgp$alt = NULL

message('SV catalog from the 1000GP ONT Vienna...')
if(!file.exists('final-vcf.unphased.vcf.gz')){
  download.file('https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/release/v1.0/final-vcf/final-vcf.unphased.vcf.gz', 'final-vcf.unphased.vcf.gz')
}
sv.ontkgpv = readSVvcf.multisamps('final-vcf.unphased.vcf.gz')
sv.ontkgpv$ref = sv.ontkgpv$alt = NULL

message('SV catalog from the 1000GP ONT Imputation...')
if(!file.exists('sniffles2_joint_sv_calls.vcf.gz')){
  download.file('http://1kgp-sv-imputation.s3-eu-west-1.amazonaws.com/sv_calls/sniffles2_joint_sv_calls.vcf.gz', 'sniffles2_joint_sv_calls.vcf.gz')
}
sv.ontkgpi = readSVvcf.multisamps('sniffles2_joint_sv_calls.vcf.gz')
sv.ontkgpi$af = sv.ontkgpi$ac / (2* (sv.ontkgpi$nrefs + sv.ontkgpi$ncalls))
sv.ontkgpi$ref = sv.ontkgpi$alt = NULL

message('SV catalog from the DGV common SVs (nstd186)...')
if(!file.exists('nstd186.GRCh38.variant_call.vcf.gz')){
  download.file('https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/nstd186.GRCh38.variant_call.vcf.gz', 'nstd186.GRCh38.variant_call.vcf.gz')
}
dgvc = readSVvcf('nstd186.GRCh38.variant_call.vcf.gz', other.field='AF')
dgvc$af = ifelse(is.na(dgvc$AF), .01, dgvc$AF)
dgvc$AF = NULL
dgvc.cnv = subset(dgvc, type=='CNV')
dgvc.cnv$type = 'DEL'
dgvc = c(dgvc, dgvc.cnv)
dgvc$type = ifelse(dgvc$type=='CNV', 'DUP', dgvc$type)
dgvc$type = ifelse(dgvc$type=="DUP", "INS", dgvc$type)
seqlevels(dgvc) = paste0('chr', seqlevels(dgvc))

message('All DGV SVs regions...')
if(!file.exists('GRCh38_hg38_variants_2020-02-25.txt')){
  download.file('http://dgv.tcag.ca/dgv/docs/GRCh38_hg38_variants_2020-02-25.txt',
                'GRCh38_hg38_variants_2020-02-25.txt')
}
dgv.heads = read.table('GRCh38_hg38_variants_2020-02-25.txt', nrows=1, header=TRUE, sep='\t')

dgva = read.table('GRCh38_hg38_variants_2020-02-25.txt', sep='\t', header=TRUE,
                  comment.char = '',
                  colClasses=c('NULL', 'character', 'integer', 'integer',
                               'NULL', 'character',
                               rep('NULL', ncol(dgv.heads)-6)))
dgva$type = ifelse(grepl('insertion', dgva$variantsubtype), 'INS', dgva$variantsubtype)
dgva$type = ifelse(grepl('deletion', dgva$type), 'DEL', dgva$type)
dgva$type = ifelse(grepl('duplication', dgva$type), 'DEL', dgva$type)
dgva$type = ifelse(dgva$type=='gain+loss', 'CNV', dgva$type)
dgva$type = ifelse(dgva$type=='gain', 'DUP', dgva$type)
dgva$type = ifelse(dgva$type=='loss', 'DEL', dgva$type)
dgva$type = ifelse(dgva$type=='inversion', 'INV', dgva$type)
dgva.cnv = subset(dgva, type=='CNV')
dgva.cnv$type = 'DEL'
dgva = rbind(dgva, dgva.cnv)
dgva$type = ifelse(dgva$type=='CNV', 'DUP', dgva$type)
dgva$type = ifelse(dgva$type=="DUP", "INS", dgva$type)
dgva = subset(dgva, type %in% c('DEL', 'INS', 'INV'))
dgva = subset(dgva, chr %in% c(1:22, 'X', 'Y'))
dgva$chr = paste0('chr', dgva$chr)
dgva$variantsubtype = NULL
dgva = makeGRangesFromDataFrame(dgva, keep.extra.columns = TRUE)
dgva$size = width(dgva)

message('Simple repeats and homopolymers...')
if(!file.exists('simpleRepeat.grch38.txt.gz')){
  download.file('https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz', 'simpleRepeat.grch38.txt.gz')
}
sr = read.table('simpleRepeat.grch38.txt.gz', as.is=TRUE, header=FALSE)
max.base.percent = apply(sr[,12:15], 1, function(x) max(x))
homopol = sr[which(max.base.percent>95),]
homopol = GRanges(homopol[,2], IRanges(homopol[, 3], homopol[, 4]))
sr = reduce(GRanges(sr[,2], IRanges(sr[, 3], sr[, 4])))
gc()

message('Segmental duplications...')
if(!file.exists('genomicSuperDups.grch38.txt.gz')){
  download.file('https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz', 'genomicSuperDups.grch38.txt.gz')
}
sd = read.table('genomicSuperDups.grch38.txt.gz', as.is=TRUE)
sd = sd[,c(2:4, 27)]
colnames(sd) = c('chr', 'start', 'end', 'fracMatch')
sd = makeGRangesFromDataFrame(sd, keep.extra.columns=TRUE)

message('Assembly gaps, centromeres, telomeres...')
if(!file.exists('gap.grch38.txt.gz')){
  download.file('https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz', 'gap.grch38.txt.gz')
}
gaps = read.table('gap.grch38.txt.gz', as.is=TRUE, sep='\t')
gaps = GRanges(gaps$V2, IRanges(gaps$V3, gaps$V4), type=gaps$V8)
if(!file.exists('centromeres.grch38.txt.gz')){
  download.file('https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz', 'centromeres.grch38.txt.gz')
}
cent = read.table('centromeres.grch38.txt.gz', as.is=TRUE, sep='\t')
gaps = c(gaps, GRanges(cent$V2, IRanges(cent$V3, cent$V4), type='centromere'))

message('OMIM...')
if(!file.exists('omim-genemap2.txt')){
  stop('Register and download omim-genemap2.txt from https://www.omim.org/downloads')
}
omim = read.table('omim-genemap2.txt', as.is=TRUE, sep='\t', comment.char='#', quote = '')
formatIn <-function(inh){
  in.pats = c(AD='autosomal dominant', AR='autosomal recessive',
           XL='x-linked', XLD='x-linked dominant', XLR='x-linked recessive',
           YL='y-linked')
  grep.o  = lapply(in.pats, function(pat){
    grepl(pat, inh, ignore.case=TRUE)
  })
  grep.o = do.call(cbind, grep.o)
  apply(grep.o, 1, function(x) paste(names(in.pats)[which(x)], collapse=';'))
}
omim = omim %>% mutate(inheritance=formatIn(V13))
## genes are listed with alternative names
omim.genes.l = strsplit(omim$V7, ', ')
omim.genes.idx = rep(1:length(omim.genes.l), unlist(lapply(omim.genes.l, length)))
omim.genes = tibble(
  gene = unlist(omim.genes.l),
  primaryName = unlist(lapply(omim.genes.l, function(x) 1:length(x) == 1)),
  inheritance=omim$inheritance[omim.genes.idx],
  omim.pheno=omim$V13[omim.genes.idx]!=''
)

message('ClinGen genes...')
if(!file.exists('Clingen-Gene-Disease-Summary.csv')){
  stop('Download Clingen-Gene-Disease-Summary CSV from https://search.clinicalgenome.org/kb/downloads')
}
clingen.genes = read.csv('Clingen-Gene-Disease-Summary.csv', skip=6, header=FALSE, as.is=TRUE)
clingen.genes = clingen.genes[,c(1,3,5,7,10)]
colnames(clingen.genes) = c("gene","disease", "inheritance","classification","gcep")
clingen.genes = unique(clingen.genes)

message('ClinGen dosage genes...')
if(!file.exists('Clingen-Dosage-Sensitivity.csv')){
  stop('Download Clingen-Dosage-Sensitivity-Summary CSV from https://search.clinicalgenome.org/kb/downloads')
}
clingen.dos.genes = read.csv('Clingen-Dosage-Sensitivity.csv', skip=6, header=FALSE, as.is=TRUE)
clingen.dos.genes = clingen.dos.genes[,c(1,3,4)]
colnames(clingen.dos.genes) = c("gene","haplo", "triplo")
haplo.ins = c("Sufficient Evidence for Haploinsufficiency",
              "Emerging Evidence for Haploinsufficiency",
              "Gene Associated with Autosomal Recessive Phenotype")
triplo.ins = c("Sufficient Evidence for Triplosensitivity",
               "Emerging Evidence for Triplosensitivity")
clingen.dos.genes = unique(clingen.dos.genes) %>% mutate(haplo=haplo %in% haplo.ins,
                                                         triplo=triplo %in% triplo.ins) %>%
  filter(haplo | triplo)

message('GenCC genes...')
if(!file.exists('gencc-submissions.csv')){
  stop('Download submissions csv from https://search.thegencc.org/download')
}
gencc.genes = read.csv('gencc-submissions.csv', header=TRUE, as.is=TRUE)
gencc.genes = gencc.genes %>%
  select(gene_symbol, disease_title, classification_title, moi_title) %>% unique

disease.genes = unique(c(clingen.genes$gene, gencc.genes$gene_symbol))

message('Genehancer...')
if(!file.exists('GeneHancer.gff.gz')){
  stop('Request access at https://www.genecards.org/Guide/Datasets')
}
df = read.table('GeneHancer.gff.gz', as.is=TRUE, comment.char='', sep='\t', header=TRUE)
colnames(df)[1] = 'seqnames'
df = df %>% select(seqnames, start, end, score, attributes)
att.l = strsplit(df$attributes, ';')
df.d = df[rep(1:nrow(df), unlist(lapply(att.l, function(x) (length(x)-1)/2))),]
rownames(df.d) = NULL
df.d$gene = gsub('connected_gene=', '', unlist(lapply(att.l, function(x) x[seq(2, length(x), 2)])))
df.d$gene_score = as.numeric(gsub('score=', '', unlist(lapply(att.l, function(x) x[seq(3, length(x), 2)]))))
df.d = df.d %>% select(-attributes)
dis.genes = unique(c(clingen.genes$gene, gencc.genes$gene_symbol))
dis.enh = df.d %>% filter(gene %in% dis.genes) %>%
  makeGRangesFromDataFrame(keep.extra.columns=TRUE)

message('ENCODE Enhancers...')
if(!file.exists('GRCh38-ELS.bed')){
  download.file('https://downloads.wenglab.org/cCREs/GRCh38-ELS.bed', 'GRCh38-ELS.bed')
}
enc.enh = read.table('GRCh38-ELS.bed', as.is=TRUE, sep='\t')
enc.enh = with(enc.enh, GRanges(V1, IRanges(V2, V3)))
prom.dis = subset(genc, type=='promoter' & gene_type=='protein_coding' & gene_name %in% disease.genes)
d.df = findOverlaps(enc.enh, prom.dis, maxgap=1e6) %>% as.data.frame %>%
  mutate(gene=prom.dis$gene_name[subjectHits],
         distance=distance(enc.enh[queryHits], prom.dis[subjectHits])) %>%
  arrange(distance) %>% group_by(queryHits) %>% 
  mutate(nearest=rep(c(TRUE, FALSE), c(1, n()-1)))
dis.enh.enc = enc.enh[d.df$queryHits,]
mcols(dis.enh.enc) = d.df[,c('gene', 'distance', 'nearest')]

message('Mappable regions...')
if(!file.exists('k100.umap.bed.gz')){
  download.file('https://bismap.hoffmanlab.org/raw/hg38/k100.umap.bed.gz', 'k100.umap.bed.gz')
}
map = read.table('k100.umap.bed.gz', as.is=TRUE, skip=1, sep='\t')
map = GRanges(map[,1], IRanges(map[,2], map[,3]))

message("JARVIS high-score regions...")
jarvis = read.table('jarvis_ranked_deduplicated.hg38.perc99.ext10.minext90.tsv.gz', as.is=TRUE, sep='\t', fill=TRUE)
jarvis = with(jarvis, GRanges(paste0('chr', V1), IRanges(V2, V3)))

message('Save in RData file...')
sv.catalogs = list(
  gnomad=gnomad,
  giraffe=sv.giraffe,
  hprc=sv.hprc,
  hgsvc=sv.hgsvc,
  DGVcommon=dgvc,
  ONTKGPhigh=sv.ontkgp,
  ONTKGPmed=sv.ontkgpv,
  ONTKGPlow=sv.ontkgpi
)

save(sv.catalogs,
     dgva,
     genc, types.ranked, types.labels,
     clinsv, dis.enh, dis.enh.enc, jarvis, 
     gene.pli, gene.phaplo, gene.ptriplo,
     omim.genes, clingen.genes, clingen.dos.genes, gencc.genes,
     sr, homopol, sd, gaps, map, file='annotation_database.RData')
