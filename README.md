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
## SV database
##### TODO: upload .RData file somewhere public

miniwdl run --as-me -i test/test.inputs.json wdl/workflow.wdl
```
