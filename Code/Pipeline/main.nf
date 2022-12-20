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
       	gunzip -c ${reference} > reference.fa
        """
}


process index_bwa {
	cpus {8}
	
	input:	
	file reference from ref_bwa_ch

	output:
	file "*" into bwa_index  

	"""
	bwa-mem2 index -p reference ${reference} 
	"""
}

process samtools_index {
	cpus {8}

	input:
	file reference from ref_sam_ch

	output:
	file("reference*.fai") into sam_known_variants, sam_extract_known, sam_filter_known_snps, sam_filter_known_indels, sam_bsqr_ch, sam_variant_ch, sam_merge_ch, sam_genotype_ch, sam_extract_ch, sam_filter_snps_ch, sam_filter_indels_ch

	"""
	samtools faidx reference.fa
	"""
}

process sequence_dictionary {
	container 'broadinstitute/gatk:latest'
	cpus {8}

	input:
	file reference from ref_dict_ch
	
	output:
	file("${reference.baseName}.dict") into dict_known_variants, dict_extract_known, dict_filter_known_snps, dict_filter_known_indels, dict_bsqr_ch, dict_variant_ch, dict_merge_ch, dict_genotype_ch, dict_extract_ch, dict_filter_snps_ch, dict_filter_indels_ch

	"""
	gatk CreateSequenceDictionary \
	-R ${reference} \
	-O ${reference.baseName}.dict
	"""
}

process align_sort {
	cpus {16}

	input:
	set sample, file(reads) from read_pair_ch  
	file "*" from bwa_index

	output:
	file ("${sample}_sorted.bam") into sorted_aligns

	"""
	bwa-mem2 mem -t 16 -R '@RG\\tID:${sample}\\tSM:${sample}\\tPL:Illumina' \
	reference ${reads.get(0)} ${reads.get(1)} | samtools sort --threads 16 -o ${sample}_sorted.bam -O BAM
        """
}


process bam_markdup {
	container 'broadinstitute/gatk:latest'
	cpus {8}
	
	input:
	file aligns from sorted_aligns

	output:	
	file("${aligns.baseName}_dedup.bam") into dedup_index_ch, dedup_known_ch, dedup_bsqr_ch, dedup_ch

	"""
	gatk MarkDuplicates \
	--INPUT ${aligns} \
	--OUTPUT ${aligns.baseName}_dedup.bam \
	-M marked_dup_metrics.txt
	"""
}

process index_bam {
	container 'staphb/samtools:latest'
	cpus {8}
	
	input:
	file reads from dedup_index_ch

	output:
	file("${reads.baseName}*.bai") into indexed_bam_ch	

	"""
	samtools index ${reads}
	"""
} 

process known_variants {
	container 'broadinstitute/gatk:latest'
	cpus {16}

	input:
	file reference from ref_known_variants_ch
	file index from sam_known_variants
	file dict from dict_known_variants
	file reads from dedup_known_ch
        file indexed_reads from indexed_bam_ch

	output:
	file("${reads.baseName}.vcf.gz") into known_variants_ch
	file("${reads.baseName}.vcf.gz.tbi") into known_index_ch

	"""
	gatk HaplotypeCaller \
	-R ${reference} \
	-O ${reads.baseName}.vcf.gz \
	-I ${reads} \
	--native-pair-hmm-threads 8 \
	--QUIET true
	"""	
}

process extract_known_snp_indel {
	container 'broadinstitute/gatk:latest'
	cpus {8}

	input:
	file reference from ref_extract_known_ch
	file index from sam_extract_known
        file dict from dict_extract_known
	file variants from known_variants_ch
	file index from known_index_ch

	output:
	file("${variants.simpleName}_snps.vcf.gz") into known_snps_ch
	file("${variants.simpleName}_indels.vcf.gz") into known_indels_ch
	file("${variants.simpleName}_snps.vcf.gz.tbi") into known_snps_index
	file("${variants.simpleName}_indels.vcf.gz.tbi") into known_indels_index
	"""
	gatk SelectVariants \
	-R ${reference} \
	-V ${variants} \
	-select-type SNP \
	-O ${variants.simpleName}_snps.vcf.gz

	gatk SelectVariants \
	-R ${reference} \
	-V ${variants} \
	-select-type INDEL \
	-O ${variants.simpleName}_indels.vcf.gz
	"""
}

process filter_known_snps {
	container 'broadinstitute/gatk:latest'
	cpus{16}
	
	input:
	file reference from ref_filter_known_snps_ch
	file index from sam_filter_known_snps
        file dict from dict_filter_known_snps
	file snps from known_snps_ch
	file snp_inddex from known_snps_index

	output:
	file("${snps.simpleName}_filtered.vcf.gz") into filtered_known_snps_ch

	"""
	gatk VariantFiltration \
        -R ${reference} \
        -V ${snps} \
	-O ${snps.simpleName}_filtered.vcf.gz \
	-filter-name "QD_filter" -filter "QD $params.qd_snp" \
	-filter-name "FS_filter" -filter "FS $params.fs_snp" \
	-filter-name "MQ_filter" -filter "MQ $params.mq_snp" \
	-filter-name "MQRankSum_filter" -filter "MQRankSum $params.mqrs_snp" \
	-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum $params.rprs_snp"
	"""
}

process filter_known_indels {
	container 'broadinstitute/gatk:latest'
	cpus{16}
	
	input:
	file reference from ref_filter_known_indels_ch
	file index from sam_filter_known_indels
        file dict from dict_filter_known_indels
	file indels from known_indels_ch
	file indels_inddex from known_indels_index

	output:
	file("${indels.simpleName}_filtered.vcf.gz") into filtered_known_indels_ch
	"""	
	gatk VariantFiltration \
	-R ${reference} \
	-V ${indels} \
	-O ${indels.simpleName}_filtered.vcf.gz \
	-filter-name "QD_filter" -filter "QD $params.qd_indel" \
	-filter-name "FS_filter" -filter "FS $params.fs_indel" \
	-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum $params.rprs_indel"	
	"""
}


process bqsr_initialization {
	container 'broadinstitute/gatk:latest'
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
	container 'broadinstitute/gatk:latest'
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
        container 'broadinstitute/gatk:latest'
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
	container 'broadinstitute/gatk:latest'
        cpus{16}

        input:
        file reference from ref_merge_ch
        file index from sam_merge_ch
        file dict from dict_merge_ch
        file variants from variants_ch
        file var_index from var_index_ch

        output:
        file("mergedGVCFs.g.vcf.gz") into merged_gvcf_ch
	file("mergedGVCFs.g.vcf.gz.tbi") into merged_index

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
        container 'broadinstitute/gatk:latest'
	cpus{16}

        input:
        file reference from ref_genotype_ch
        file index from sam_genotype_ch
        file dict from dict_genotype_ch
        file merged_var from merged_gvcf_ch
	file variant_index from merged_index

        output:
        file {variants.vcf.gz} into complete_vcf_ch
	file {variants.vcf.gz.tbi} into complete_vcf_index

        """
         gatk --java-options "-Xmx4g" GenotypeGVCFs \
        -R ${reference} \
        -V ${merged_var} \
        -O variants.vcf.gz
        """
}

process extract_snp_indel {
	container 'broadinstitute/gatk:latest'
	cpus {16}

	input:
	file reference from ref_extract_ch
	file index from sam_extract_ch
	file dict from dict_extract_ch
	file variants from complete_vcf_ch
	file variant_index from complete_vcf_index

	output:
	file("snps_recal.vcf.gz") into snps_ch
	file("indels_recal.vcf.gz") into indels_ch
	file("snps_recal.vcf.gz.tbi") into snps_recal_index
        file("indels_recal.vcf.gz.tbi") into indels_recal_index

	"""
	gatk SelectVariants \
	-R ${reference} \
	-V ${variants} \
	-select-type SNP \
	-O snps_recal.vcf.gz
	gatk SelectVariants \
	-R ${reference} \
	-V ${variants} \
	-select-type INDEL \
	-O indels_recal.vcf.gz
	"""
}
process filter_snps {
	container 'broadinstitute/gatk:latest'
	cpus{16}
	
	publishDir "${outdir}", mode: 'copy'

	input:
	file reference from ref_filter_snps_ch
        file index from sam_filter_snps_ch
        file dict from dict_filter_snps_ch
	file snps from snps_ch
	file snp_index from snps_recal_index

	output:
	file("snps/filtered_snps.vcf.gz") into snp_for_snpeff
	
	"""
	mkdir snps
	gatk VariantFiltration \
        -R ${reference} \
        -V ${snps} \
	-O snps/filtered_snps.vcf.gz \
	-filter-name "QD_filter" -filter "QD $params.qd_snp" \
	-filter-name "FS_filter" -filter "FS $params.fs_snp" \
	-filter-name "MQ_filter" -filter "MQ $params.mq_snp" \
	-filter-name "MQRankSum_filter" -filter "MQRankSum $params.mqrs_snp" \
	-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum $params.rprs_snp"
	"""
}
process filter_indels {
	container 'broadinstitute/gatk:latest'
	cpus{16}
	
	input:
	file reference from ref_filter_indels_ch
	file index from sam_filter_indels_ch
        file dict from dict_filter_indels_ch
	file indels from indels_ch
	file indels_index from indels_recal_index

	output:
	file("filtered_indels.vcf.gz") into indels_for_snpeff

	"""
	gatk VariantFiltration \
	-R ${reference} \
	-V ${indels} \
	-O filtered_indels.vcf.gz \
	-filter-name "QD_filter" -filter "QD $params.qd_indel" \
	-filter-name "FS_filter" -filter "FS $params.fs_indel" \
	-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum $params.rprs_indel"
	"""
}

process annotate_snp{
	container 'dceoy/snpeff:latest'
	cpus{8}
	
	publishDir "${outdir}", mode: 'copy', pattern: '*.html'
	publishDir "${outdir}", mode: 'copy', pattern: '*.ann.vcf'
	input:
	file filtered_snp from snp_for_snpeff

	output:
	file ("ann_snp.ann.vcfz")
	file("snp_annotation_report.html")

	"""
	java -Xmx8g -jar snpEff.jar -v -stats snp_annotation_report.html Sorghum_bicolor ${filtered_snp}  > ann_snp.ann.vcfz
	"""
}

process annotate_indels{
	container 'dceoy/snpeff:latest'
        cpus{8}

        publishDir "${outdir}", mode: 'copy', pattern: '*.html'
        publishDir "${outdir}", mode: 'copy', pattern: '*.ann.vcf'
        input:
        file filtered_indels from indels_for_snpeff

        output:
        file ("ann_indels.ann.vcfz")
        file("indel_annotation_report.html")

        """
        java -Xmx8g -jar snpEff.jar -v -stats indel_annotation_report.html Sorghum_bicolor ${filtered_snp}  > ann_indels.ann.vcfz
        """
}

