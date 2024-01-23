version 1.0

workflow denovo_svs {

    meta {
        author: "Jean Monlong"
        email: "jmonlong@ucsc.edu"
        description: "Annotate candidate de novo SVs"
    }

    parameter_meta {
        VCF: "Input VCF from the child. Can be gzipped/bgzipped."
        VCF_P1: "Input VCF from parent 1. Can be gzipped/bgzipped."
        VCF_P2: "Input VCF from parent 2. Can be gzipped/bgzipped."
        SV_DB_RDATA: "RData file with the databases used for SV annotation (e.g. SV catalogs, simple repeat annotation)."
        BAM: "Sorted and indexed BAM with the long reads for the child. Optional. If present, SVs will be re-genotyped to help filtering out false-positives."
        BAI: "Index for BAM of the child"
        BAM_P1: "Sorted and indexed BAM with the long reads for parent 1. Optional. If present, SVs will be re-genotyped to help filtering out false-positives."
        BAI_P1: "Index for BAM of parent 1"
        BAM_P2: "Sorted and indexed BAM with the long reads for parent 1. Optional. If present, SVs will be re-genotyped to help filtering out false-positives."
        BAI_P2: "Index for BAM of parent 1"
        REFERENCE_FASTA: "Reference fasta. Optional. Used for SV in silico validation when BAM and BAI are provided."
    }
    
    input {
        File VCF
        File VCF_P1
        File VCF_P2
        File SV_DB_RDATA
        File? BAM
        File? BAI
        File? BAM_P1
        File? BAI_P1
        File? BAM_P2
        File? BAI_P2
        File? REFERENCE_FASTA
    }

    # annotate de novo candidate SVs (in child but not parents) with allele frequency
    call annotate_denovo_svs {
        input:
        input_vcf=VCF,
        input_vcf_p1=VCF_P1,
        input_vcf_p2=VCF_P2,
        sv_db_rdata=SV_DB_RDATA
    }
    
    # regenotype SVs with local pangenomes using vg to provide some in silico "validation"
    if (defined(BAM) && defined(BAI) && defined(BAM_P1) && defined(BAI_P1) && defined(BAM_P2) && defined(BAI_P2) && defined(REFERENCE_FASTA)){
        call genotype_svs_trio_with_vg {
            input:
            input_vcf=annotate_denovo_svs.vcf,
            bam=BAM,
            bai=BAI,
            bam_p1=BAM_P1,
            bai_p1=BAI_P1,
            bam_p2=BAM_P2,
            bai_p2=BAI_P2,
            reference_fasta=REFERENCE_FASTA
        }
    }

    File final_vcf = select_first([genotype_svs_trio_with_vg.vcf_stringent, annotate_denovo_svs.vcf])
        
    output {
        File vcf = final_vcf
        File? vcf_lenient = genotype_svs_trio_with_vg.vcf_lenient
    }
}

task annotate_denovo_svs {
    input {
        File input_vcf
        File input_vcf_p1
        File input_vcf_p2
        File sv_db_rdata
        Int memSizeGB = 8
        Int threadCount = 2
        Int diskSizeGB = 5*round(size(input_vcf, "GB") + size(sv_db_rdata, 'GB')) + 30
    }

    String basen = sub(sub(basename(input_vcf), ".vcf.bgz$", ""), ".vcf.gz$", "")
    
    command <<<
        set -eux -o pipefail

        # extract SVs for trio
        bcftools view -i "STRLEN(REF)>=40 | MAX(STRLEN(ALT))>=40" -Oz -o svs.vcf.gz ~{input_vcf}
        bcftools view -i "STRLEN(REF)>=40 | MAX(STRLEN(ALT))>=40" -Oz -o svs.p1.vcf.gz ~{input_vcf_p1}
        bcftools view -i "STRLEN(REF)>=40 | MAX(STRLEN(ALT))>=40" -Oz -o svs.p2.vcf.gz ~{input_vcf_p2}

        # annotate SVs
        Rscript /opt/scripts/annotate_denovo_svs.R svs.vcf.gz svs.p1.vcf.gz svs.p2.vcf.gz ~{sv_db_rdata} svs.denovo.vcf
        
        # sort and gzip VCF
        bcftools sort -Oz -o ~{basen}.svs_denovo.vcf.gz svs.denovo.vcf
    >>>

    output {
        File vcf = "~{basen}.svs_denovo.vcf.gz"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/jmonlong/svannotate_sveval@sha256:c694bd3db95a07a49a231856c4662326976cc7f1d9875f746adb122ffc094ad6"
        preemptible: 1
    }
}

task genotype_svs_trio_with_vg {
    input {
        File input_vcf
        File? bam
        File? bai
        File? bam_p1
        File? bai_p1
        File? bam_p2
        File? bai_p2
        File? reference_fasta
        Int memSizeGB = 8
        Int threadCount = 16
        Int diskSizeGB = 5*round(size(input_vcf, "GB") + size(bam, 'GB') + size(bam_p1, 'GB') + size(bam_p2, 'GB') + size(reference_fasta, 'GB')) + 30
    }

    String basen = sub(sub(basename(input_vcf), ".vcf.bgz$", ""), ".vcf.gz$", "")
    
    command <<<
        set -eux -o pipefail

        ## link and index reference fasta
        ln -s ~{reference_fasta} ref.fa
        samtools faidx ref.fa

        ## link BAM files and indexes to make sure the index is found
        ln -s ~{bam} reads.bam
        ln -s ~{bai} reads.bam.bai
        ln -s ~{bam_p1} reads.p1.bam
        ln -s ~{bai_p1} reads.p1.bam.bai
        ln -s ~{bam_p2} reads.p2.bam
        ln -s ~{bai_p2} reads.p2.bam.bai

        ## genotype only rare SVs (common SVs are unlikely to be real de novo)
        bcftools view -i "AF<0.01" -o input.vcf.gz -Oz ~{input_vcf}
        
        ## run the validation script in parallel across smCores cores
        python3 /opt/scripts/validate-svs-trio.py -b reads.bam,reads.p1.bam,reads.p2.bam -f ref.fa -v ~{input_vcf} -t ~{threadCount} -o ~{basen}.svval.vcf.gz

        ## extract two sets of de novo candidates
        ## stringent: some evidence in child, no evidence in parents
        bcftools view -i "RS_PROP_P1==0 && RS_PROP_P2==0 && RS_PROP>0" -o ~{basen}.denovo.stringent.vcf.gz -Oz output.vcf.gz
        ## lenient: some evidence in child, low evidence in parents
        bcftools view -i "RS_PROP_P1<0.1 && RS_PROP_P2<0.1 && RS_PROP>0" -o ~{basen}.denovo.lenient.vcf.gz -Oz output.vcf.gz
    >>>

    output {
        File vcf_stringent = "~{basen}.denovo.stringent.vcf.gz"
        File vcf_lenient = "~{basen}.denovo.lenient.vcf.gz"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/jmonlong/svvalidate_vgcall@sha256:599f73d3471acdecab0e9107488f21a55f35b18805dc71fe937961ad2a1157df"
        preemptible: 1
    }
}
