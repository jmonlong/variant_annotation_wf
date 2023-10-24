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
        SV_DB_RDATA: "RData file with the databases used for SV annotation (e.g. SV catalogs, dbVar clinical SVs, DGV)."
        SORT_INDEX_VCF: "Should the output VCF be sorted, bgzipped, and indexed? Default: true"
        BAM: "Sorted and indexed BAM with the long reads for the child. Optional. If present, SVs will be re-genotyped to help filtering false-positives out."
        BAI: "Index for BAM of the child"
        BAM_P1: "Sorted and indexed BAM with the long reads for parent 1. Optional. If present, SVs will be re-genotyped to help filtering false-positives out."
        BAI_P1: "Index for BAM of parent 1"
        BAM_P2: "Sorted and indexed BAM with the long reads for parent 1. Optional. If present, SVs will be re-genotyped to help filtering false-positives out."
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
        Boolean SORT_INDEX_VCF = true
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

    File sv_vcf = select_first([genotype_svs_trio_with_vg.vcf, annotate_denovo_svs.vcf])

    # sort annotated VCF
    if (SORT_INDEX_VCF){
        call sort_vcf {
            input:
            input_vcf=sv_vcf
        }
    }
    
    File final_vcf = select_first([sort_vcf.vcf, sv_vcf])
    
    output {
        File vcf = final_vcf
        File? vcf_index = sort_vcf.vcf_index
    }
}

task sort_vcf {
    input {
        File input_vcf
        Int memSizeGB = 4
        Int diskSizeGB = 5*round(size(input_vcf, "GB")) + 20
    }

    String basen = sub(sub(basename(input_vcf), ".vcf.bgz$", ""), ".vcf.gz$", "")
    
    command <<<
    set -eux -o pipefail

    bcftools sort -Oz -o ~{basen}.formatted.vcf.gz ~{input_vcf}
    bcftools index -t -o ~{basen}.formatted.vcf.gz.tbi ~{basen}.formatted.vcf.gz
    >>>
    
    output {
        File vcf = "~{basen}.formatted.vcf.gz"
        File? vcf_index = "~{basen}.formatted.vcf.gz.tbi"
    }
    
    runtime {
        memory: memSizeGB + " GB"
        cpu: 1
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/biocontainers/bcftools@sha256:f3a74a67de12dc22094e299fbb3bcd172eb81cc6d3e25f4b13762e8f9a9e80aa"
        preemptible: 1
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

        # extract SVs and small variants
        bcftools view -i "STRLEN(REF)>=30 | MAX(STRLEN(ALT))>=30" -Oz -o svs.vcf.gz ~{input_vcf}
        bcftools view -i "STRLEN(REF)>=30 | MAX(STRLEN(ALT))>=30" -Oz -o svs.p1.vcf.gz ~{input_vcf_p1}
        bcftools view -i "STRLEN(REF)>=30 | MAX(STRLEN(ALT))>=30" -Oz -o svs.p2.vcf.gz ~{input_vcf_p2}

        # annotate SVs
        Rscript /opt/scripts/annotate_denovo_svs.R svs.vcf.gz svs.p1.vcf.gz svs.p2.vcf.gz ~{sv_db_rdata} svs.denovo.vcf
        
        # merge back SVs
        bcftools sort -Oz -o ~{basen}.svs_denovo.vcf.gz svs.denovo.vcf
    >>>

    output {
        File vcf = "~{basen}.svs_denovo.vcf.gz"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/jmonlong/svannotate_sveval@sha256:c4a0ac3dd9176cfc0bec81980d254fa18eb62e0d390d2d01b54c57441a86a65c"
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
    Int smCores = threadCount / 2
    
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

        ## run the validation script in parallel across smCores cores
        ln -s ~{input_vcf} input.vcf.gz
        snakemake --snakefile /opt/scripts/Snakefile_trio --cores ~{smCores}

        mv output.vcf.gz ~{basen}.svval.vcf.gz
    >>>

    output {
        File vcf = "~{basen}.svval.vcf.gz"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/jmonlong/svvalidate_vgcall@sha256:15b9b2837c57f9d44452b73ea94b8cce54387eb65d52dd3f22251c09dc1a66ba"
        preemptible: 1
    }
}
