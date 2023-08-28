Container to validate SVs by re-genotyping a SV candidate with vg

The [`validate-svs.py`](validate-svs.py) script is included (available at `/opt/scripts/validate-svs.py` in the container).
It will construct a local pangenome for each SV, map long reads with minigraph, and genotype them with `vg pack/call`.

The inputs of `validate-svs.py` are:

```
root@d012132f9326:/home# python3 /opt/scripts/validate-svs.py --help
usage: validate-svs.py [-h] -b B -f F -v V [-d D] [-o O] [-t T]

optional arguments:
  -h, --help  show this help message and exit
  -b B        BAM file (indexed)
  -f F        reference FASTA file (indexed)
  -v V        variants in VCF (can be bgzipped)
  -d D        output directory
  -o O        output (annotated) VCF (will be bgzipped if ending in .gz)
  -t T        number of threads used by the tools (vg and minigraph))
```

Two new INFO fields are added in the output VCF:

- `RS_PROP` with the proportion of supporting reads.
- `RS_AD` with the read support for the reference and alternate alleles (e.g. `3,5` for 3 reference-supporting reads and 5 SV-supporting reads).

Latest version: `quay.io/jmonlong/svvalidate_vgcall:0.5` 

To build locally and upload to [quay.io](https://quay.io/jmonlong/svvalidate_vgcall).
