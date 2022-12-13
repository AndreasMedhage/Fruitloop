// Paths
input_reads = params.inputreads
ref_genome = params.refgenome
out_dir = params.outdir


read_pair_ch = Channel.fromFilePairs( "${input_reads}/$params.pattern", type: 'file')

ref_ch = Channel.fromPath("${ref_genome}").collect()

// unzip reference
process unzip_reference {
        cpus {1}

        input:
        file reference from ref_ch

        output:
        file("reference.fa") into ref_bwa_ch, ref_sam_ch, ref_dict_ch, ref_known_variants_ch, ref_extract_known_ch, ref_filter_known_snps_ch, ref_filter_known_indels_ch, ref_bsqr_ch, ref_variant_ch, ref_merge_ch, ref_genotype_ch, ref_extract_ch, ref_filter_snps_ch, ref_filter_indels_ch

        """
	if (file ${reference} | grep -q compressed ); then
        	gunzip -c ${reference} > reference.fa
	else
		cp ${reference} reference.fa
	fi
        """
}


process index_bwa {
	cpus {8}
	
	input:	
	file ref from ref_bwa_ch

	output:
	file "*" into bwa_index  

	"""
	bwa-mem2 index -p reference ${ref} 
	"""
}

process samtools_index {
	cpus {8}

	input:
	file ref from ref_sam_ch

	output:
	file("reference*.fai") into sam_known_variants, sam_extract_known, sam_filter_known_snps, sam_filter_known_indels, sam_bsqr_ch, sam_varaiant_ch, sam_merge_ch, sam_genotype_ch, sam_extract_ch, sam_filter_snps_ch, sam_filter_indels_ch

	"""
	samtools faidx reference.fa
	"""
}

process sequence_dictionary {
	cpus {8}

	input:
	file ref from ref_dict_ch
	
	output:
	file("${ref.baseName}.dict") into dict_known_variants, dict_extract_known, dict_filter_known_snps, dict_filter_known_indels, dict_bsqr_ch, dict_variant_ch, dict_merge_ch, dict_genotype_ch, dict_extract_ch, dict_filter_snps_ch, dict_filter_indels_ch

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
	file("${aligns.baseName}_dedup.bam") into dedup_index_ch, dedup_known_ch, dedup_bsqr_ch, dedup_ch

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
	file reads from dedup_index_ch

	output:
	file("${reads}*.bai") into indexed_bam_ch	

	"""
	samtools index ${reads}
	"""
} 

process known_variants {
	cpus {16}

	input:
	file reference from ref_known_variants_ch
	file index from sam_known_variants
	file dict from dict_known_variants
	file reads from dedup__known_ch
        file indexed_reads from indexed_bam_ch

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

process extract_known_snp_indel {
	cpus {8}

	input:
	file reference from ref_extract_known_ch
	file index from sam_extract_known
        file dict from dict_extract_known
	file variants from known_variants_ch

	output:
	file("*_snps.vcf") into known_snps_ch
	file("*_indels.vcf") into known_indels_ch
	"""
	gatk SelectVariants \
	-R ${reference} \
	-V ${variants} \
	-select-type SNP \
	-O ${variants.baseName}_snps.vcf

	gatk SelectVariants \
	-R ${reference} \
	-V ${variants} \
	-select-type INDEL \
	-O ${variants.baseName}_indels.vcf
	"""
}

process filter_known_snps {
	cpus{16}
	
	input:
	file reference from ref_filter_known_snps_ch
	file index from sam_filter_known_snps
        file dict from dict_filter_known_snps
	file snps from known_snps_ch

	output:
	file("*_filtered.vcf") into filtered_known_snps_ch

	"""
	gatk VariantFiltration \
        -R ${reference} \
        -V ${snps} \
	-O ${snps.baseName}_filtered.vcf \
	-filter-name "QD_filter" -filter "QD $params.qd_snp" \
	-filter-name "FS_filter" -filter "FS $params.fs_snp" \
	-filter-name "MQ_filter" -filter "MQ $params.mq_snp" \
	-filter-name "MQRankSum_filter" -filter "MQRankSum $params.mqrs_snp" \
	-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum $params.rprs_snp"
	"""
}

process filter_known_indels {
	cpus{16}
	
	input:
	file reference from ref_filter_known_indels_ch
	file index from sam_filter_known_indels
        file dict from dict_filter_known_indels
	file indels from known_indels_ch

	output:
	file("*_filtered.vcf") into filtered_known_indels_ch
	"""	
	gatk VariantFiltration \
	-R ${reference} \
	-V ${indels} \
	-O ${indels.baseName}_filtered.vcf \
	-filter-name "QD_filter" -filter "QD $params.qd_indel" \
	-filter-name "FS_filter" -filter "FS $params.fs_indel" \
	-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum $params.rprs_indel"	
	"""
}


process bqsr_initialization {
        cpus {16}

        input:
        file reference from ref_bsqr_ch
        file index from sam_bsqr_ch
        file dict from dict_bsqr_ch
	file marked_dedup from dedup_bsqr_ch
        file known_snps from filtered_known_snps_ch
        file known_indels from filtered_known_indels_ch

        output:
        file("${marked_dedup.baseName}_recal_data.table") into baserecal_table_ch
        file("${marked_dedup.baseName}_recal_data.table") into baserecal_report_ch

        """
	gatk IndexFeatureFile \
	-I ${known_snps}
	
	gatk IndexFeatureFile \
        -I ${known_indels}

        gatk BaseRecalibrator \
        -I ${marked_dedup} \
        -O ${marked_dedup.baseName}_recal_data.table \
        -R ${reference} \
        --known-sites ${known_snps} \
        --known-sites ${known_indels}
        """
}

process apply_bqsr {
        cpus {16}

        input:
        file bam_dedup from dedup_ch
        file table from baserecal_table_ch

        output:
        file("${bam_dedup.baseName}_bqsr.bam") into bam_bqsr_ch
	file("${bam_dedup.baseName}_bqsr.bai") into bam_bqsr_index_ch
        """
        gatk ApplyBQSR \
        -I ${bam_dedup} \
        -bqsr ${bam_dedup.baseName}_recal_data.table \
        -O ${bam_dedup.baseName}_bqsr.bam \
        -OBI true
        """
}

process variant_calling {
        cpus {16}

        publishDir "${outdir}", mode: 'copy', pattern: '*.g.vcf.gz'
        publishDir "${outdir}", mode: 'copy', pattern: '*_final.bam'

        input:
        file reference from ref_variant_ch
        file index from sam_variant_ch
        file dict from dict_variant_ch
	file reads from bam_bqsr_ch
        file indexed_reads from bam_bqsr_index_ch


        output:
        file("GVCF/${reads.baseName}.g.vcf.gz") into variants_ch
        file("GVCF/${reads.baseName}.g.vcf.gz.tbi") into var_index_ch
        file("final_bam/${reads.baseName}_final.bam") into final_bam_ch
		
	"""
        mkdir GVCF
        mkdir final_bam

        gatk HaplotypeCaller \
        -R ${reference} \
        -I ${reads} \
        -O GVCF/${reads.baseName}.g.vcf.gz \
        --native-pair-hmm-threads 8 \
        -ERC GVCF \
        -bamout final_bam/${reads.baseName}_final.bam \
        --QUIET true
        """
}

process merge_gvcf {
        cpus{16}

        input:
        file reference from ref_merge_ch
        file index from sam_merge_ch
        file dict from dict_merge_ch
        file variants from variants_ch
        file var_index from var_index_ch

        output:
        file("mergedGVCFs.g.vcf.gz") into merged_gvcf_ch

        """
        for vcf in \$(ls *.vcf.gz); do
                echo \$vcf >> variant_files.list
        done

        gatk CombineGVCFs \
        -R ${reference} \
        --variant variant_files.list \
        -O mergedGVCFs.g.vcf.gz
        """
}

process genotype_gvcf {
        cpus{16}

        input:
        file reference from ref_genotype_ch
        file index from sam_genotype_ch
        file dict from dict_genotype_ch
        file merged_var from merged_gvcf_ch

        output:
        file {variants.vcf.gz}

        """
         gatk --java-options "-Xmx4g" GenotypeGVCFs \
        -R ${reference} \
        -V ${merged_var} \
        -O variants.vcf.gz
        """
}

process extract_snp_indel {
	cpus {16}

	input:
	file reference from ref_extract_ch
	file index from sam_extract_ch
	file dict from dict_extract_ch
	file variants from merged_gvcf_ch

	output:
	file("snps_recal.vcf") into snps_ch
	file("indels_recal.vcf") into indels_ch

	"""
	gatk SelectVariants \
	-R ${reference} \
	-V ${variants} \
	-select-type SNP \
	-O snps_recal.vcf
	gatk SelectVariants \
	-R ${reference} \
	-V ${variants} \
	-select-type INDEL \
	-O indels_recal.vcf
	"""
}
process filter_snps {
	cpus{16}
	
	publishDir "${outdir}", mode: 'copy'

	input:
	file reference from ref_filter_snps_ch
        file index from sam_filter_snps_ch
        file dict from dict_filter_snps_ch
	file snps from snps_ch

	output:
	file("snps/filtered_snps.vcf") into filtered_snps.ch

	"""
	mkdir snps
	gatk VariantFiltration \
        -R ${reference} \
        -V ${snps} \
	-O snps/filtered_snps.vcf \
	-filter-name "QD_filter" -filter "QD $params.qd_snp" \
	-filter-name "FS_filter" -filter "FS $params.fs_snp" \
	-filter-name "MQ_filter" -filter "MQ $params.mq_snp" \
	-filter-name "MQRankSum_filter" -filter "MQRankSum $params.mqrs_snp" \
	-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum $params.rprs_snp"
	"""
}
process filter_indels {
	cpus{16}
	
	input:
	file reference from ref_filter_indels_ch
	file index from sam_filter_indels_ch
        file dict from dict_filter_indels_ch
	file indels from indels_ch

	output:
	file("filtered_indels.vcf") into filtered_indels.ch

	"""
	gatk VariantFiltration \
	-R ${reference} \
	-V ${indels} \
	-O filtered_indels.vcf \
	-filter-name "QD_filter" -filter "QD $params.qd_indel" \
	-filter-name "FS_filter" -filter "FS $params.fs_indel" \
	-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum $params.rprs_indel"
	"""
}
