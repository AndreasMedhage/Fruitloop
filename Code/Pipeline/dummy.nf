QD_snp=params.qd-snp
FS_snp=params.fs-snp
MQ_snp=params.mq-snp
SOR_snp=params.sor-snp
MQRS_snp=params.mqrs-snp
RPRS_snp=params.rprs-snp
QD_indel=params.qd-indel
FS_indel=params.fs-idenl
SOR_indel=params.sor-indel

process variant_calling {
        cpus {16}

	publishDir "${outdir}", mode: 'copy', pattern: '*.g.vcf.gz'
	publishDir "${outdir}", mode: 'copy', pattern: '*_final.bam' 

        input:
        file reads from bam_bqsr_ch
        file indexed_reads from indexed_bam_ch
        file reference from unzip_ref_5_ch
        file index from sam_index_3
        file dict from seq_dict_3

        output:
        file("GVCF/${reads.baseName}.g.vcf.gz") into variants_ch
        file("${reads.baseName}*.idx*") into var_index_ch
	file("final_bam/*_final.bam") into final_bam_ch

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
	file reference from unzip_ref_6_ch
        file index from sam_index_4
        file dict from seq_dict_4
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
	file reference from unzip_ref_7_ch
        file index from sam_index_5
        file dict from seq_dict_5
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
