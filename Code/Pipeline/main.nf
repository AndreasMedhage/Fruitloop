trimmed_reads = params.trimmedreads
ref_genome = params.refgenome
out_dir = params.outdir
readpair_pattern = params.pattern
known_variants = params.knownvar

read_pair_ch = Channel.fromFilePairs( "${trimmed_reads}/${readpair_pattern}", type: 'file')

ref_ch = Channel.fromPath("${ref_genome}").collect()


process index_bwa {
	cpus {8}
	
	input:	
	file ref from ref_ch

	output:
	file "*" into bwa_index  

	"""
	bwa-mem2 index -p reference ${ref} 
	"""
}

process samtools_index {
	cpus {8}

	input:
	file ref from ref_ch

	output:
	file("reference*.fai") into sam_index

	"""
	zcat ${ref} | bgzip -c > reference.bgz
	samtools faidx reference.bgz
	"""
}

process sequence_dictionary {
	cpus {8}

	input:
	file ref from ref_ch
	
	output:
	file("${ref}.dict") into seq_dict

	"""
	java -jar $PICARD_ROOT/picard.jar \
	CreateSequenceDictionary \
	R=${ref} \
	O=${ref}.dict
	"""
}

process align_bwa {
	cpus {8}


	input:
	set sample, file(reads) from read_pair_ch  
	file "*" from bwa_index

	output:
	file "${sample}.sam" into align_ch

	"""
	bwa-mem2 mem -t 8 -R '@RG\\tID:${sample}\\tSM:${sample}\\tPL:Illumina' \
	reference ${reads.get(0)} ${reads.get(1)} > ${sample}.sam
        """
}

process bam_sort {
	cpus {8}

	input: 	
	file aligns from align_ch
	
	output:
	file ("${aligns.baseName}_sorted.bam") into sorted_aligns

	"""
	samtools sort -o ${aligns.baseName}_sorted.bam -O BAM ${aligns}
	"""
}

process bam_markdup {
	cpus {8}
	
	input:
	file aligns from sorted_aligns

	output:	
	file("${aligns.baseName}_dedup_reads.bam") into dedup_ch, bam_qc_ch

	"""
	java -jar $PICARD_ROOT/picard.jar MarkDuplicates \
	--INPUT ${aligns} \
	--OUTPUT ${aligns.baseName}_dedup_reads.bam \
	-M marked_dup_metrics.txt
	"""
}

