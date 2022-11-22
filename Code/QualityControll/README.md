QualityControl - directory storing code used for trimming and qualty controll

	initialize.sh - bash script for running pipeline
		indir= path directory containging forward and reverse reads
		outdir1= path outdirectory trimming
		outdir2= path outdirectory quality control 
		pattern_rev= pattern distinguish forward and reverse reads ex "*R[1,2]*" for the Uppsala sequences

	main.nf - main nextflow file for the for quality controll pipeline
		--rawreads - is the input directory, containing forward and reverse reads marked with "R1" and "R2" respectivly
                --outdirtrimmed - is the outdirectory for fast qc. Will create folders "trimmed" containing trimmed sequences and "report_files" containing fast p quality reports as json files
                --outdirqc - is the outdirectory for quality controll report
		--readpairpattern - pattern to distinguish forward and reverse read
		-c - config file
		-profile uppmax - nextflow profile
                --project <project code>
                --clusterOptions "-M <server>"


	nextflow.config - configuration file for main.nf
		process.executor ='slurm'
		process.time="<time> h"
		process.cpus="<cpus>"
		process.memory="<memory> GB"
		process.clusterOptions - allows further parameters

To run:

change indir, outdir1 and outdir2 in initialize.sh  
run the script "sh initialize.sh"

