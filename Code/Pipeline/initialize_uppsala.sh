#!/bin/bash -l

module load bioinfo-tools Nextflow/20.10.0 picard/2.27.5 bwa-mem2/2.2.1-20211213-edc703f samtools/1.14 GATK/4.3.0.0 bcftools/1.14

export NXF_HOME="/proj/snic2022-23-502/private/Fruitloop/Code/Pipeline"

nextflow run main.nf --trimmedreads '/proj/snic2022-23-502/private/Fruitloop/Data/Trimmed/testuppsala/trimmed_files'\
 --refgenome "/proj/snic2022-23-502/private/Fruitloop/Data/Sequences/Reference/GCF*.fna.gz"\
 --outdir "/proj/snic2022-23-502/private/Fruitloop/Res/Pipeline/test/"\
 --pattern "*R[1,2]*.gz"\
 -c nextflow.config\
 -profile uppmax --project "snic2022-22-1048" --clusterOptions "-M rackham" -resume


