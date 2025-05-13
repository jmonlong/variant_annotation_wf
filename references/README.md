## Misc annotations in R object

The [`prepare-annotations.R`](prepare-annotations.R) script downloads the publicly accessible files or will stop to tell you to download the files that require registration (e.g. GeneHancer, OMIM), or run some commands.

Or you can prepare all the data necessary in advance by following the instructions below.

Download publicly available files:

```sh
## Gene annotations
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gtf.gz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz
wget https://zenodo.org/record/6347673/files/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz?download=1
## "Giraffe", HPRC and HGSVC SV catalogs with frequency annotation
wget https://cgl.gi.ucsc.edu/data/giraffe/calling/vggiraffe-sv-2504kgp-svsites.tsv.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.0-mc.grch38.vcfbub.a100k.wave.vcf.gz
bcftools view -i "STRLEN(REF)>30 | MAX(STRLEN(ALT))>30" hprc-v1.0-mc.grch38.vcfbub.a100k.wave.vcf.gz | gzip > hprc-v1.0-mc.grch38.vcfbub.a100k.wave.svs.vcf.gz
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/release/v1.0/integrated_callset/freeze3.sv.alt.vcf.gz
## ONT 1000GP
wget https://s3.amazonaws.com/1000g-ont/Gustafson_etal_2024_preprint_SUPPLEMENTAL/20240423_jasmine_intrasample_noBND_custom_suppvec_alphanumeric_header_JASMINE.vcf.gz
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/release/v1.0/final-vcf/final-vcf.unphased.vcf.gz
wget http://1kgp-sv-imputation.s3-eu-west-1.amazonaws.com/sv_calls/sniffles2_joint_sv_calls.vcf.gz
## DGV/dbVar catalogs
wget https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/nstd186.GRCh38.variant_call.vcf.gz
wget http://dgv.tcag.ca/dgv/docs/GRCh38_hg38_variants_2020-02-25.txt
wget https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/nstd166.GRCh38.variant_call.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/tsv/nstd102.GRCh38.variant_call.tsv.gz
## UCSC Genome Browser annotations
wget -O simpleRepeat.grch38.txt.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz
wget -O genomicSuperDups.grch38.txt.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz
wget -O gap.grch38.txt.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz
## ENCODE
wget https://downloads.wenglab.org/cCREs/GRCh38-ELS.bed
## mappability
wget https://bismap.hoffmanlab.org/raw/hg38/k100.umap.bed.gz
```

Some files require registration or are updated regularly (so you might want to control/log when they are downloaded):

- Register and download `omim-genemap2.txt` from https://www.omim.org/downloads
- Download Clingen *Gene-Disease Validity* and *Dosage Sensitivity Curations (Genes Only)* CSVs from https://search.clinicalgenome.org/kb/downloads
    - Rename (or link) to `Clingen-Gene-Disease-Summary.csv` and `Clingen-Dosage-Sensitivity.csv`.
- Download `gencc-submissions.csv` *Submissions (csv)* from https://search.thegencc.org/download
- Request access to GeneHancer at https://www.genecards.org/Guide/Datasets
    - Rename (or link) the GFF file to `GeneHancer.gff.gz`.


For the [JARVIS scores](https://www.nature.com/articles/s41467-021-21790-4), we first downloaded the files from [https://az.app.box.com/v/jarvis-gwrvis-scores/folder/159706242087](https://az.app.box.com/v/jarvis-gwrvis-scores/folder/159706242087), and unzipped the `JARVIS.zip`.
Then, we then use a naive sliding-window approach to find regions with high JARVIS scores, implemented in the `find-high-score-regions.py`.

```sh
for CHR in `seq 1 22`
do
    unzip -p JARVIS.zip JARVIS/jarvis_ranked_deduplicated.${CHR}.hg38.tsv.gz | zcat | python3 find-high-score-regions.py
done | gzip > jarvis_ranked_deduplicated.hg38.perc99.ext10.minext90.tsv.gz
```

After downloading and preparing all those files, we can run the R script that compiles everything in one R object:

```sh
Rscript prepare-annotations.R
```

If you would like to skip this and directly get the file that we use for GRCh38, contact us.

## Tandem repeat expansion in public controls

We compute the length of annotated repeats by combining SV calls for each haplotypes using the [`prepare_tr_from_vcf.py`](prepare_tr_from_vcf.py).
We compare the length of those repeat expansions/contractions to the distribution in control samples, here HPRC and HGSVC2 samples. 
This catalog of comparable repeat variation was computed by running the same script on the calls from HPRC/HGSVC and produce `hprc-hgsvc-tr-catalog.tsv.gz`.

```sh
## download the HPRC/HGSVC calls if not done already above
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.0-mc.grch38.vcfbub.a100k.wave.vcf.gz
bcftools view -i "STRLEN(REF)>30 | MAX(STRLEN(ALT))>30" hprc-v1.0-mc.grch38.vcfbub.a100k.wave.vcf.gz | gzip > hprc-v1.0-mc.grch38.vcfbub.a100k.wave.svs.vcf.gz
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/release/v1.0/integrated_callset/freeze3.sv.alt.vcf.gz

## prepare the repeat length variation files
zcat hprc-v1.0-mc.grch38.vcfbub.a100k.wave.svs.vcf.gz | python3 ../scripts/prepare_tr_from_vcf.py | gzip > hprc-v1.0-mc.grch38.vcfbub.a100k.wave.svs.sum.tsv.gz
zcat hgsvc.freeze3.sv.alt.vcf.gz | python3 ../scripts/prepare_tr_from_vcf.py | gzip > hgsvc.freeze3.sv.alt.sum.tsv.gz

## combine
Rscript prepare-tr-catalog.R
gzip hprc-hgsvc-tr-catalog.tsv
```
