#!/bin/bash -l
#SBATCH -A snic2022-22-1048
#SBATCH -M rackham
#SBATCH -p node
#SBATCH -n 2
##SBATCH -t 25:00:00
#SBATCH -J fastp_uppsala_seq
#SBATCH --mail-type=ALL
#SBATCH --mail-user kaavya.venkateswaran.5938@student.uu.se

#to go to the folder where the script is present
cd /crex/proj/snic2022-23-502/private/Fruitloop/Code/QualityControll

#load required modules and set the environment
module load bioinfo-tools Nextflow/20.10.0 fastp/0.23.2 MultiQC/1.12

export NXF_HOME="/proj/snic2022-23-502/private/Fruitloop/Code/QualityControll"

#run the nextflow
nextflow run main.nf --rawreads "/proj/snic2022-23-502/private/Fruitloop/Data/Sequences/Uppsala/*.gz" \
--outdirtrimmed "/proj/snic2022-23-502/private/Fruitloop/Res/QualityControll/Uppsala" \
--outdirqc "/proj/snic2022-23-502/private/Fruitloop/Data/Trimmed/Uppsala" -c nextflow.config \
-profile uppmax --project "snic2022-22-1048" --clusterOptions "-M rackham" -resume


