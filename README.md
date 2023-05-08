# Variant annotation

Simple variant annotation.

- Variants are annotated with their predicted impact on the genes using SNPeff
- Then with public variant databases:
    - Small variants are annotated with with gnomAD's allele frequencies, and ClinVar's clinical significance.
    - Structural variants are annotated with the frequency in SV catalogs, and their overlap with dbVar clinical SVs, or DGV.

## Test locally

```
## download GRCh38.105 database
wget https://snpeff.blob.core.windows.net/databases/v5_1/snpEff_v5_1_GRCh38.105.zip

miniwdl run --as-me -i test/test.inputs.json wdl/workflow.wdl
```

## Prepare annotation files

### gnomAD

To reduce the size of the VCF with the variants in gnomAD, keeping only the field that we want to use (AF):

```sh
for CHR in `seq 1 22` X Y
do
    curl https://storage.googleapis.com/gcp-public-data--gnomad/release/3.0/vcf/genomes/gnomad.genomes.r3.0.sites.chr${CHR}.vcf.bgz | zcat | bcftools annotate -x ^INFO/AF,ID,FILTER,QUAL -e 'FILTER~"AC0"' -Oz -o gnomad_chr.genomes.r3.0.sites.chr${CHR}.vcf.gz
done
bcftools concat gnomad_chr.genomes.r3.0.sites.chr*.vcf.gz | bcftools sort -Oz -o gnomad.genomes.r3.0.sites.small.vcf.bgz -m 3G
rm gnomad_chr.genomes.r3.0.sites.chr*.vcf.gz
```

The final bgzipped VCF is only 4.5G. 

### dbNSFP

To reduce the size of the dbNSFP database used by SNPeff:

```sh
## download dbNSFP from the SNPeff repo
wget -O dbNSFP4.3.zip https://usf.box.com/shared/static/9r6iamhldji4c3vodwebh3947vgrvsng
## Uncompress
unzip dbNSFP4.3.zip
## Create a single file version
zcat dbNSFP4.3a_variant.chr1.gz | head -n 1 > header.txt
zcat dbNSFP4.3a_variant.chr* | grep -v "^#" >> header.txt
bgzip header.txt > dbNSFP4.3.txt.gz
## (optional) list columns in full txt file. helps finding the column numbers 
zcat dbNSFP4.3.txt.gz | head -1 | awk 'BEGIN{RS="\t"}{N=N+1;print N" "$0}' | less
## keep only some columns
zcat dbNSFP4.3.txt.gz | cut -f 1-4,77-79,118-120,155-157 | bgzip > dbNSFP4.3.small.txt.gz
## index with tabix
tabix -b 2 -e 2 -s 1 dbNSFP4.3.small.txt.gz
```

If we keep only information about GERP, CADD and MetaRNN, the file size decreases from 35G to about 2G.
