trimmed_reads = params.trimmedreads
ref_genome = params.refgenome
out_dir = params.outdir
readpair_pattern = params.pattern

read_pair_ch = Channel.fromFilePairs( "${trimmed_reads}/${readpair_pattern}", type: 'file')

ref_ch = Channel.fromPath("${ref_genome}").collect()


process unzip_reference {
        cpus {1}

        input:
        file reference from ref_ch

        output:
        file("reference.fa") into unzip_ref_ch, unzip_ref_2_ch, unzip_ref_3_ch, unzip_ref_4_ch, unzip_ref_5_ch

        """
        gunzip -c ${reference} > reference.fa
        """
}


process index_bwa {
	cpus {8}
	
	input:	
	file ref from unzip_ref_ch

	output:
	file "*" into bwa_index  

	"""
	bwa-mem2 index -p reference ${ref} 
	"""
}

process samtools_index {
	cpus {8}

	input:
	file ref from unzip_ref_2_ch

	output:
	file("reference*.fai") into sam_index, sam_index_2

	"""
	samtools faidx reference.fa
	"""
}

process sequence_dictionary {
	cpus {8}

	input:
	file ref from unzip_ref_3_ch
	
	output:
	file("${ref.baseName}.dict") into seq_dict, seq_dict_2

	"""
	java -jar $PICARD_ROOT/picard.jar \
	CreateSequenceDictionary \
	-R ${ref} \
	-O ${ref.baseName}.dict
	"""
}

process align_bwa {
	cpus {16}

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
	file("${aligns.baseName}_dedup.bam") into dedup_ch, dedup_2_ch, dedup_3_ch, dedup_4_ch

	"""
	java -jar $PICARD_ROOT/picard.jar MarkDuplicates \
	--INPUT ${aligns} \
	--OUTPUT ${aligns.baseName}_dedup.bam \
	-M marked_dup_metrics.txt
	"""
}

process index_bam {
	cpus {8}
	
	input:
	file reads from dedup_ch

	output:
	file("${reads}*.bai") into indexed_bam_ch	

	"""
	samtools index ${reads}
	"""
} 

process known_variants {
	cpus {16}

	input:
	file reads from dedup_2_ch
	file indexed_reads from indexed_bam_ch
	file reference from unzip_ref_4_ch
	file index from sam_index
	file dict from seq_dict

	output:
	file("${reads.baseName}.vcf") into known_variants_ch
	file("${reads.baseName}.vcf.idx") into index_ch

	"""
	gatk HaplotypeCaller \
	-R ${reference} \
	-O ${reads.baseName}.vcf \
	-I ${reads} \
	--native-pair-hmm-threads 8 \
	--QUIET true
	"""	
}

process merge_known_variants {
	cpus {8}

	input:
	file ('*.vcf') from known_variants_ch.collect()
	file ('*.vcf.idx') from index_ch.collect()

	output:
	file("known_variants.vcf") into known_vcf_ch 
	file("known_variants.vcf.idx") into vcf_index_ch	

	"""
	for vcf in \$(ls *vcf); do
		echo \$vcf >> variant_files.list
	done
	
	bcftools --file-list variant_files.list --output known_variants.vcf --threads 8	
	"""
}

process bqsr_initialization {
        cpus {8}

        publishDir "${out_dir}", pattern: '*.table'

        input:
        file marked_dedup from dedup_3_ch
        file reference from unzip_ref_5_ch
        file index from sam_index_2
        file dict from seq_dict_2
        file known_vcf from known_vcf_ch
        file known_index from vcf_index_ch

        output:
        file("bqs_report/${marked_dedup.baseName}_recal_data.table") into baserecal_table_ch
        file("bqs_report/data.table") into baserecal_report_ch

        """
        mkdir bqs_report
        gatk BaseRecalibrator \
        -I ${marked_dedup} \
        -O bqs_report/${marked_dedup.baseName}_recal_data.table \
        -R ${reference} \
        --known-sites ${known_vcf}
        """
}

process apply_bqsr {
        cpus {8}

        publishDir "${out_dir}", mode: 'copy', overwrite: false

        input:
        file bam_dedup from dedup_4_ch
        file table from baserecal_table_ch

        output:
        file("bqsr/${marked_dedup.baseName}_bqsr.bam") into bam_bqsr_ch

        """
        mkdir bqsr
        gatk ApplyBQSR \
        -I ${bam_dedup} \
        -bqsr ${bam_dedup.baseName}_recal_data.table \
        -O bqsr/${marked_dedup.baseName}_bqsr.bam
        """
}

