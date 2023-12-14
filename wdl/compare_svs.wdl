version 1.0

workflow compare_svs {

    meta {
        author: "Jean Monlong"
        email: "jmonlong@ucsc.edu"
        description: "Compare SVs in two VCFs, e.g. from different callers"
    }

    parameter_meta {
        VCF_1: "Input VCF 1. Can be gzipped/bgzipped."
        VCF_2: "Input VCF 2. Can be gzipped/bgzipped."
        SV_DB_RDATA: "RData file with the databases used for SV annotation (e.g. SV catalogs, simple repeat annotation)."
    }
    
    input {
        File VCF_1
        File SV_DB_RDATA
    }

    # annotate de novo candidate SVs (in child but not parents) with allele frequency
    call compare_svs_vcf {
        input:
        input_vcf_1=VCF_1,
        input_vcf_2=VCF_2,
        sv_db_rdata=SV_DB_RDATA
    }
    
    File final_vcf = select_first([genotype_svs_trio_with_vg.vcf_stringent, annotate_denovo_svs.vcf])
        
    output {
        File ann_vcf_1 = compare_svs_vcf.vcf_1
        File ann_vcf_2 = compare_svs_vcf.vcf_2
    }
}

task compare_svs_vcf {
    input {
        File input_vcf_1
        File input_vcf_2
        File sv_db_rdata
        String sample_name_1 = ''
        String sample_name_2 = ''
        Int memSizeGB = 8
        Int threadCount = 2
        Int diskSizeGB = 5*round(size(input_vcf_1, "GB") + size(sv_db_rdata, 'GB')) + 30
    }

    String basen1 = sub(sub(basename(input_vcf_1), ".vcf.bgz$", ""), ".vcf.gz$", "")
    String basen2 = sub(sub(basename(input_vcf_2), ".vcf.bgz$", ""), ".vcf.gz$", "")
    
    command <<<
        set -eux -o pipefail

        # extract SVs for trio
        bcftools view -i "STRLEN(REF)>=40 | MAX(STRLEN(ALT))>=40" -Oz -o svs.1.vcf.gz ~{input_vcf_1}
        bcftools view -i "STRLEN(REF)>=40 | MAX(STRLEN(ALT))>=40" -Oz -o svs.2.vcf.gz ~{input_vcf_2}

        # annotate SVs
        Rscript /opt/scripts/annotate_sv_overlap.R svs.1.vcf.gz svs.2.vcf.gz ~{sv_db_rdata} svs.1.ann.vcf svs.2.ann.vcf ~{sample_name_1} ~{sample_name_1}
        
        # sort and gzip VCF
        bcftools sort -Oz -o ~{basen1}.comp.vcf.gz svs.1.ann.vcf
        bcftools sort -Oz -o ~{basen2}.comp.vcf.gz svs.2.ann.vcf
    >>>

    output {
        File ann_vcf_1 = "~{basen1}.comp.vcf.gz"
        File ann_vcf_2 = "~{basen2}.comp.vcf.gz"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/jmonlong/svannotate_sveval@sha256:c694bd3db95a07a49a231856c4662326976cc7f1d9875f746adb122ffc094ad6"
        preemptible: 1
    }
}
