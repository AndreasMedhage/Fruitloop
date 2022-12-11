#!/bin/bash -l

module load bioinfo-tools Nextflow/20.10.0 picard/2.27.5 bwa-mem2/2.2.1-20211213-edc703f samtools/1.14 GATK/4.3.0.0

export NXF_HOME="/proj/snic2022-23-502/private/Fruitloop/Code/Pipeline"

nextflow run main.nf --trimmedreads '/proj/snic2022-23-502/private/Fruitloop/Data/Trimmed/testuppsala/trimmed_files'\
 --refgenome "/proj/snic2022-23-502/private/Fruitloop/Data/Sequences/Reference/GCF*.fna.gz"\
 --outdir "/proj/snic2022-23-502/private/Fruitloop/Res/Pipeline/test/"\
 --pattern "*R[1,2]*.gz"\
 --qd_snp "< 2.0" --fs_snp "> 60.0 " --mq_snp "< 40.0" --mqrs_snp "< -12.5" --rprs_snp "< 8.0" --qd_indel "< 2.0" --fs_indel "> 200.0 " --rprs_indel "> 20.0"\
 -c nextflow.config\
 -profile uppmax --project "snic2022-22-1048" --clusterOptions "-M rackham" -resume



