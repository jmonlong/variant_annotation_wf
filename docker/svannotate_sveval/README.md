Container for the SV annotation task. 
It has [bcftools](https://github.com/samtools/bcftools) (v1.17), and R with the [sveval](https://github.com/jmonlong/sveval) package.

The [`annotate_svs.R`](annotate_svs.R) annotation script is also included (available at `/opt/scripts/annotate_svs.R` in the container).

Latest version: `quay.io/jmonlong/svannotate_sveval:0.1` 

To build locally and upload to [quay.io](https://quay.io/repository/jmonlong/svannotate_sveval):

```
docker build -t jmonlong-svannotate_sveval:0.1 .
docker tag jmonlong-svannotate_sveval:0.1 quay.io/jmonlong/svannotate_sveval:0.1
docker push quay.io/jmonlong/svannotate_sveval:0.1
```
