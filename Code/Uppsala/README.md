Pipline - directory storing code of the main pipeline

        initialize.sh - bash script for running pipeline
                

        main.nf - main nextflow file for the for quality controll pipeline
		--inputreads - input directory, containing forward and reverse reads
		--refgenome - reference genome to map reads to
		--outdir - out directory, where output files should be pubished
		-c - config file
 		-profile - nextflow profile
		--project - project code
		--clusterOptions - computing cluster options, for ex  "-M <server name>"
		-resume - To resume last run"

        nextflow.config - configuration file for main.nf
                process.executor ='slurm'
                process.time="<time> h"
                process.cpus="<cpus>"
                process.memory="<memory> GB"
                process.clusterOptions - allows further parameters

		params {
			pattern - pattern for destinguising forward and revers, for example "*R[1,2]*.gz" forward marked with R1 and reverse marked
			qd_snp -  QualByDepth filter for SNPs, ex "< 2.0"
			fs_snp - FisherStrand (Fisher’s Exact Test) filter for SNPs, ex "> 60.0 "
			mq_snp = Maping Quality Root Mean Square filter for SNPs, ex "< 40.0"
			mqrs_snp = Mann-Whitney Rank Sum Test (mapping quality) filter for SNPs, ex "< -12.5"
			rprs_snp = Mann-Whitney Rank Sum Test (distance end of read) filter for SNPs, ex "< 8.0"
			qd_indel = QualByDepth filter for indels, ex "< 2.0"
			fs_indel = FisherStrand (Fisher’s Exact Test) filter for indels, ex "> 200.0 "
			rprs_indel = Mann-Whitney Rank Sum Test (distance end of read) filter for indels, ex "> 20.0"
		}

To run:

edit inizialize.sh
	set trimmedreads to path to direcory containging forward and reverse reads
	set outdirectory to path directory to output into
	set refgenome to path to reference genome file

edit config file set up pattern forward/reverse and filter parameter values

run the script "sh initialize.sh"


