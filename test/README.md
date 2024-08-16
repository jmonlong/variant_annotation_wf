# Test dataset

Small test data used only to validate the syntax/scripts of the workflow.

It's a subset of the human genome `chr7:97357275-104567465` with simulated reads.
Run the following *required* section before running the tests.

## Required: download large files

Some files are trickier to downsize so we still need to download them.
That's the case for:

- the reference genome FASTA
- SnpEff database
- AnnotSV database

They are not provided in this repo. 
See the code chunk below to prepare them.

```sh
# download reference genome FASTA:
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
samtools faidx hg38.fa

# SnpEff
wget https://snpeff.blob.core.windows.net/databases/v5_1/snpEff_v5_1_GRCh38.105.zip
unzip snpEff_v5_1_GRCh38.105.zip

# AnnotSV
docker run -it -v `pwd`:/app -w /app -u `id -u $USER` quay.io/jmonlong/annotsv:3.4
# - # IN DOCKER CONTAINER
cp -r /build/AnnotSV .
cd AnnotSV
PREFIX=. make install-human-annotation
mv share/AnnotSV ../AnnotSV_annotations
cd ..
rm -fr AnnotSV
# - # OUT DOCKER CONTAINER
```

## Provided: downsized files

The rest are provided in this repo and are tiny files relevant for just the `chr7:97357275-104567465` slice of the genome.
For info, this data was prepared with:

```sh
# subset some annotation files (annotation_database.RData, hprc-hgsvc-tr-catalog.tsv.gz)
Rscript prepare_test.R

# simulare some reads
python3 prepare_test.py -f hg38.fa -r chr7:97357200-104567500 -v test.vcf.gz -n 1000 -l 10000
# map and index them
minimap2 hg38.fa reads.fastq -a | samtools sort -o reads.bam
samtools index reads.bam
```

# Run tests

## Snakemake workflow

Assuming you're in the root of the repository, run the workflow with:

```sh
snakemake --use-singularity --configfile test/snakemake.config.yaml -p reports --cores 4
```

The config file points to all the relevant files (FASTA, databases, sample information).
It also points to a TSV file with sample information, i.e. the path to the variants (VCF) and reads (BAM) for each sample.
See [snakemake.config.yaml](snakemake.config.yaml) and [snakemake.sample_info.tsv](snakemake.sample_info.tsv) for examples of those two files.

## WDL workflow

*Soon...*
