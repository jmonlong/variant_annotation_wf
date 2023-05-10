version 1.0

workflow annotate_variants {

    meta {
        author: "Jean Monlong"
        email: "jmonlong@ucsc.edu"
        description: "Annotate a VCF with SNPeff, frequencies in gnomAD, and ClinVar"
    }

    parameter_meta {
        VCF: "Input VCF. Can be gzipped/bgzipped."
        SNPEFF_DB: "SNPeff annotation bundle (zip file)"
        SNPEFF_DB_NAME: "Name of SNPeff annotation bundle. E.g. GRCh38.105"
        GNOMAD_VCF: "VCF with all variants in gnomAD and their allele frequency (AF INFO field). Must be sorted, bgzipped, and indexed."
        GNOMAD_VCF_INDEX: "Index for GNOMAD_VCF (.tbi file)."
        CLINVAR_VCF: "VCF with all variants in ClinVar and their clinical significance (CLNSIG INFO field). Must be sorted, bgzipped, and indexed."
        CLINVAR_VCF_INDEX: "Index for CLINVAR_VCF (.tbi file)."
        SV_DB_RDATA: "RData file with the databases used for SV annotation (e.g. SV catalogs, dbVar clinical SVs, DGV)."
        SPLIT_MULTIAL: "Should multiallelic variants be split into biallelic records? Default: true"
        SORT_INDEX_VCF: "Should the output VCF be sorted, bgzipped, and indexed? Default: true"
        DBNSFP_DB: "dbNSFP annotation bundle (block-gzip)"
        DBNSFP_DB_INDEX: "Index for DBNSFP_DB"
        BAM: "Sorted and indexed BAM with the long reads. Optional. If present, SVs will be re-genotyped to help filtering false-positives out."
        BAM_INDEX: "Index for BAM"
        REFERENCE_FASTA: "Reference fasta. Optional. Used for SV in silico validation when BAM and BAM_INDEX are provided."
    }
    
    input {
        File VCF
        File? SNPEFF_DB
        String? SNPEFF_DB_NAME
        File? GNOMAD_VCF
        File? GNOMAD_VCF_INDEX
        File? CLINVAR_VCF
        File? CLINVAR_VCF_INDEX
        File? SV_DB_RDATA
        File? DBNSFP_DB
        File? DBNSFP_DB_INDEX
        File? BAM
        File? BAM_INDEX
        File? REFERENCE_FASTA
        Boolean SPLIT_MULTIAL = true
        Boolean SORT_INDEX_VCF = true
    }

    # split multi-allelic variants
    if (SPLIT_MULTIAL){
        call split_multiallelic_vcf {
            input:
            input_vcf=VCF
        }
    }

    # annotate variants with predicted effect based on gene annotation
    File current_vcf = select_first([split_multiallelic_vcf.vcf, VCF])

    if(defined(SNPEFF_DB) && defined(SNPEFF_DB)){
        call annotate_with_snpeff {
            input:
            input_vcf=current_vcf,
            snpeff_db=select_first([SNPEFF_DB]),
            db_name=select_first([SNPEFF_DB_NAME])
        }
    }

    File annotated_vcf = select_first([annotate_with_snpeff.vcf, current_vcf])

    # annotate SNVs/indels with frequency in gnomAD and presence in ClinVar
    # note: first filter variants to keep those with high/moderate impact or with predicted loss of function (speeds up DB matching a lot)
    if (defined(GNOMAD_VCF) && defined(GNOMAD_VCF_INDEX) && defined(CLINVAR_VCF) && defined(CLINVAR_VCF_INDEX) && defined(DBNSFP_DB) && defined(DBNSFP_DB_INDEX)){
        call subset_annotate_smallvars_with_db {
            input:
            input_vcf=annotated_vcf,
            gnomad_vcf=select_first([GNOMAD_VCF]),
            gnomad_vcf_index=select_first([GNOMAD_VCF_INDEX]),
            clinvar_vcf=select_first([CLINVAR_VCF]),
            clinvar_vcf_index=select_first([CLINVAR_VCF_INDEX]),
            dbnsfp_db = select_first([DBNSFP_DB]),
            dbnsfp_db_index = select_first([DBNSFP_DB_INDEX])
        }
    }
    
    File small_annotated_vcf = select_first([subset_annotate_smallvars_with_db.vcf, annotated_vcf])
    
    # annotate SVs with frequency in SV databases and presence in dbVar Clinical SVs
    if (defined(SV_DB_RDATA)){
        call annotate_sv_with_db {
            input:
            input_vcf=small_annotated_vcf,
            sv_db_rdata=select_first([SV_DB_RDATA])
        }
    }
    
    File sv_annotated_vcf = select_first([annotate_sv_with_db.vcf, small_annotated_vcf])

    # regenotype SVs with local pangenomes using vg to provide some in silico "validation"
    if (defined(BAM) && defined(BAM_INDEX) && defined(REFERENCE_FASTA)){
        call validate_svs_with_vg {
            input:
            input_vcf=sv_annotated_vcf,
            bam=select_first([BAM]),
            bam_index=select_first([BAM_INDEX]),
            reference_fasta=select_first([REFERENCE_FASTA])
        }
    }

    File sv_validated_vcf = select_first([validate_svs_with_vg.vcf, sv_annotated_vcf])

    # sort annotated VCF
    if (SORT_INDEX_VCF){
        call sort_vcf {
            input:
            input_vcf=sv_validated_vcf
        }
    }
    
    File final_vcf = select_first([sort_vcf.vcf, sv_validated_vcf])
    
    output {
        File vcf = final_vcf
        File? vcf_index = sort_vcf.vcf_index
    }
}

task split_multiallelic_vcf {
    input {
        File input_vcf
        Int memSizeGB = 4
        Int threadCount = 1
        Int diskSizeGB = 5*round(size(input_vcf, "GB")) + 20
    }

    String basen = sub(sub(basename(input_vcf), ".vcf.bgz$", ""), ".vcf.gz$", "")
    
    command <<<
    set -eux -o pipefail
    
    bcftools norm -m -both --threads ~{threadCount} -Oz -o ~{basen}.norm.vcf.gz ~{input_vcf}
    >>>
    
    output {
        File vcf = "~{basen}.norm.vcf.gz"
    }
    
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/biocontainers/bcftools@sha256:f3a74a67de12dc22094e299fbb3bcd172eb81cc6d3e25f4b13762e8f9a9e80aa"
        preemptible: 1
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

task annotate_with_snpeff {
    input {
        File input_vcf
        File snpeff_db
        String db_name
        Int memSizeGB = 16
        Int threadCount = 2
        Int diskSizeGB = 5*round(size(input_vcf, "GB") + size(snpeff_db, 'GB')) + 20
    }

    Int snpeffMem = if memSizeGB < 6 then 2 else memSizeGB - 4
    String basen = sub(sub(basename(input_vcf), ".vcf.bgz$", ""), ".vcf.gz$", "")
    
    command <<<
        set -eux -o pipefail
        
        unzip ~{snpeff_db}
        
        zcat ~{input_vcf} | snpEff -Xmx~{snpeffMem}g -nodownload -no-intergenic \
               -dataDir "${PWD}/data" ~{db_name} | gzip > ~{basen}.snpeff.vcf.gz
           >>>
           
           output {
               File vcf = "~{basen}.snpeff.vcf.gz"
           }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/biocontainers/snpeff@sha256:7ac091da707f5d63f307eef4ee57c3f0e94eed49f86bbdace3d4be3a514ed410"
        preemptible: 1
    }
}

task subset_annotate_smallvars_with_db {
    input {
        File input_vcf
        File gnomad_vcf
        File gnomad_vcf_index
        File clinvar_vcf
        File clinvar_vcf_index
        File dbnsfp_db
        File dbnsfp_db_index
        Int memSizeGB = 16
        Int threadCount = 2
        Int diskSizeGB = 5*round(size(input_vcf, "GB") + size(gnomad_vcf, 'GB') + size(clinvar_vcf, 'GB') + size(dbnsfp_db, 'GB')) + 30
    }

    Int snpsiftMem = if memSizeGB < 6 then 2 else memSizeGB - 4
    String basen = sub(sub(basename(input_vcf), ".vcf.bgz$", ""), ".vcf.gz$", "")
    
    command <<<
        set -eux -o pipefail

        ## link the database VCF to make sure their indexes can be found
        ln -s ~{gnomad_vcf} gnomad.vcf.bgz
        ln -s ~{gnomad_vcf_index} gnomad.vcf.bgz.tbi
        ln -s ~{clinvar_vcf} clinvar.vcf.bgz
        ln -s ~{clinvar_vcf_index} clinvar.vcf.bgz.tbi
        ln -s ~{dbnsfp_db} dbnsfp.txt.gz
        ln -s ~{dbnsfp_db_index} dbnsfp.txt.gz.tbi

        ## filter variants to keep those with high/moderate impact or with predicted loss of function
        ## then annotate with their frequency in gnomAD
        zcat ~{input_vcf} | SnpSift -Xmx1g filter "(ANN[*].IMPACT has 'HIGH') | (ANN[*].IMPACT has 'MODERATE') | ((exists LOF[*].PERC) & (LOF[*].PERC > 0.9))" | \
            SnpSift -Xmx~{snpsiftMem}g annotate -noId -v gnomad.vcf.bgz | gzip > ~{basen}.gnomad.vcf.gz

        ## annotate IDs with clinvar IDs and add the CLNSIG INFO field
        zcat ~{basen}.gnomad.vcf.gz | SnpSift -Xmx~{snpsiftMem}g annotate -info CLNSIG -v clinvar.vcf.bgz | gzip > ~{basen}.gnomad.clinvar.vcf.gz
        
        ## annotate IDs with dbNSFP prediction scores and conservation scores
        zcat ~{basen}.gnomad.clinvar.vcf.gz > ~{basen}.gnomad.clinvar.vcf
        SnpSift -Xmx~{snpsiftMem}g dbnsfp -v -db dbnsfp.txt.gz -f GERP++_RS,CADD_raw,CADD_phred,MetaRNN_score,MetaRNN_pred ~{basen}.gnomad.clinvar.vcf | gzip > ~{basen}.gnomad.clinvar.dbnsfp.vcf.gz
    >>>

    output {
        File vcf = "~{basen}.gnomad.clinvar.dbnsfp.vcf.gz"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/biocontainers/snpsift@sha256:049babfac841d15a92d8febfc10a25f5aa109c9fe6670af35ea79583a1c78402"
        preemptible: 1
    }
}

task annotate_sv_with_db {
    input {
        File input_vcf
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
        bcftools view -i "STRLEN(REF)<30 & MAX(STRLEN(ALT))<30" ~{input_vcf} | bcftools sort -Oz -o smallvars.vcf.gz

        # annotate SVs
        Rscript /opt/scripts/annotate_svs.R svs.vcf.gz ~{sv_db_rdata} svs.annotated.vcf
        
        # merge back SVs
        bcftools sort -Oz -o svs.annotated.vcf.gz svs.annotated.vcf
        bcftools index -t svs.annotated.vcf.gz
        bcftools index -t smallvars.vcf.gz
        bcftools concat -a -Oz -o ~{basen}.svannotated.vcf.gz smallvars.vcf.gz svs.annotated.vcf.gz
    >>>

    output {
        File vcf = "~{basen}.svannotated.vcf.gz"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/jmonlong/svannotate_sveval@sha256:ca67367a1848938c05ce7508b8a9f2a9bc7fd8c9d33c06f634075677b7b529be"
        preemptible: 1
    }
}

task validate_svs_with_vg {
    input {
        File input_vcf
        File bam
        File bam_index
        File reference_fasta
        Int memSizeGB = 8
        Int threadCount = 4
        Int diskSizeGB = 5*round(size(input_vcf, "GB") + size(bam, 'GB') + size(reference_fasta, 'GB')) + 30
    }

    String basen = sub(sub(basename(input_vcf), ".vcf.bgz$", ""), ".vcf.gz$", "")
    
    command <<<
        set -eux -o pipefail

        ## link and index reference fasta
        ln -s ~{reference_fasta} ref.fa
        samtools faidx ref.fa

        ## link BAM file and index to make sure the index is found
        ln -s ~{bam} reads.bam
        ln -s ~{bam_index} reads.bam.bai
        
        # annotate SVs
        mkdir temp
        python3 /opt/scripts/validate-svs.py -b reads.bam -f ref.fa -v ~{input_vcf} -d temp -o ~{basen}.svval.vcf.gz -t ~{threadCount}
    >>>

    output {
        File vcf = "~{basen}.svval.vcf.gz"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/jmonlong/svvalidate_vgcall@sha256:f05abdf22765837084a3f046b6fb8971c8bbc50d9f60da659dfccb7fcbaf2da7"
        preemptible: 1
    }
}
