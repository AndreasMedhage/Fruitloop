Pipline - directory pipline code

        initialize.sh - bash script for running pipeline

        main.nf - main nextflow file for the for quality controll pipeline
		--inputreads - input directory, containing forward and reverse reads
		--refgenome - reference genome to map reads to
		--outdir - out directory, where output files should be pubished
		--annotation_db - database for annotation
		-c - config file
 		-profile - nextflow profile
		--project - project code
		--clusterOptions - computing cluster options, for ex  "-M <server name>"
		-resume - To resume last run

		REWUIRE: nextflow, samtools, BWA-mem2

        nextflow.config - configuration file for main.nf
                process.executor - executor
                process.time - time allocated to each process
                process.clusterOptions - allows further parameters
		process.scratch - use scratch dir to minimize storage of intermediate files
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

	main_untested.nf - untested version of main including attempt att including gvcf files from previous run
		--inputreads - input directory, containing forward and reverse reads
                --refgenome - reference genome to map reads to
                --outdir - out directory, where output files should be pubished
		--previousgvcf - input directory containing gvcf files from previous run for inclution in join genotyping
                --annotation_db - database for annotation
                -c - config file
                -profile - nextflow profile
                --project - project code
                --clusterOptions - computing cluster options, for ex  "-M <server name>"
                -resume - To resume last run
To run:


1. edit inizialize.sh
	1.1. set trimmedreads to path to direcory containging forward and reverse reads
	1.2. set outdirectory to path directory to output into
	1.3. set refgenome to path to reference genome file
	1.4. set annotation_db

2. edit config file set up pattern
	2.1. select pattern forward/reverse
	2.2. set up filter parameter values (parameters based on standar parameters of best practice, if others parameters are needed they need to be set up in the main file)

3. module load tmux

4. Create new session
	tmux session -s "name of session"

5. Detach from session
	CTRL^B + D

6. Attach to session
	tmux attatch -t "name of session"

7. run the script "sh initialize.sh"


