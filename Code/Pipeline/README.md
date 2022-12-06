Pipline - directory storing code of the main pipeline

        initialize.sh - bash script for running pipeline
                

        main.nf - main nextflow file for the for quality controll pipeline
                --rawreads - is the input directory, containing forward and reverse reads marked with "R1" and "R2" respectivly
		--reference
                
		--trimmedreads
		--refgenome
		--outdir
		--pattern



        nextflow.config - configuration file for main.nf
                process.executor ='slurm'
                process.time="<time> h"
                process.cpus="<cpus>"
                process.memory="<memory> GB"
                process.clusterOptions - allows further parameters

To run:

change indir, outdir1 and outdir2 in initialize.sh
run the script "sh initialize.sh"


