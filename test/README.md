Small test data used only to validate the syntax/scripts of the workflow.

This data was prepared with:

```sh
Rscript prepare_test.R

python3 prepare_test.py -f ../../../../genomedata/hg38.fa -r chr7:97357200-104567500 -v test.vcf.gz -n 1000 -l 10000

minimap2 ../../../../genomedata/hg38.fa reads.fastq -a | samtools sort -o reads.bam

samtools index reads.bam
```
