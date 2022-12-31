#!/bin/bash -l


module load bioinfo-tools Nextflow/20.10.0 bwa-mem2/2.2.1-20211213-edc703f samtools/1.14 snpEff/5.1 #GATK/4.3.0.0

export NXF_HOME="/proj/snic2022-23-502/private/Fruitloop/Code/Pipeline"

nextflow run main.nf --inputreads '/proj/snic2022-23-502/private/Fruitloop/Data/Trimmed/test'\
 --refgenome "/proj/snic2022-23-502/private/Fruitloop/Data/Sequences/Reference/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna.gz"\
 --annotation_db "Sorghum_bicolor"\
 --outdir "/proj/snic2022-23-502/private/Fruitloop/Res/Pipeline/test/"\
 -c nextflow.config -profile uppmax --project "snic2022-22-1048" --clusterOptions "-M rackham"


