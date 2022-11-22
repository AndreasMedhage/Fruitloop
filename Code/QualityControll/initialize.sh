#!/bin/bash -l

module load bioinfo-tools Nextflow/20.10.0 fastp/0.23.2 MultiQC/1.12 tmux

indir="/proj/snic2022-23-502/private/Fruitloop/Data/Sequences/Novogen"
outdir1="/proj/snic2022-23-502/private/Fruitloop/Data/Trimmed/Novogen"
outdir2="/proj/snic2022-23-502/private/Fruitloop/Res/QualityControll/Novogen"
pattern_rev="*_[1,2].fastq.gz"

export NXF_HOME="/proj/snic2022-23-502/private/Fruitloop/Code/QualityControll"

nextflow run main.nf --rawreads $indir --outdirtrimmed  $outdir1 --outdirqc $outdir2 --readpairpattern $pattern_rev \
-c nextflow.config -profile uppmax --project "snic2022-22-1048" --clusterOptions "-M rackham" -resume

