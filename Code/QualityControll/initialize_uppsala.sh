#!/bin/bash -l

module load bioinfo-tools Nextflow/20.10.0 fastp/0.23.2 MultiQC/1.12

export NXF_HOME="/proj/snic2022-23-502/private/Fruitloop/Code/QualityControll"

nextflow run main.nf --rawreads '/proj/snic2022-23-502/private/Fruitloop/Data/Sequences/Uppsala' \
--outdirqc "/proj/snic2022-23-502/private/Fruitloop/Res/QualityControll/Uppsala" \
--outdirtrimmed "/proj/snic2022-23-502/private/Fruitloop/Data/Trimmed/Uppsala" -c nextflow.config \
-profile uppmax --project "snic2022-22-1048" --clusterOptions "-M rackham" 

