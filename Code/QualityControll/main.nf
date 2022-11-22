#!/usr/bin/env nextflow

raw_reads = params.rawreads
out_dir_trimmed = params.outdirtrimmed
out_dir_qc = params.outdirqc
readpair_pattern = params.readpairpattern

read_pair = Channel.fromFilePairs( "${raw_reads}/${readpair_pattern}", type: 'file')


process runFastP {
        cpus { 8 }
        publishDir "${out_dir_trimmed}/", mode: 'copy', overwrite: false
        input:
                set sample, file(in_fastq) from read_pair

        output:
                file("trimmed_files/${sample}_trimmed_*.gz") into trimmed_channel
                file("report_files/${sample}*.json") into json_report
        """
        mkdir trimmed_files
        mkdir report_files
        fastp \
        -i ${in_fastq.get(0)} -I ${in_fastq.get(1)} \
        -o trimmed_files/${sample}_trimmed_R1.fq.gz -O trimmed_files/${sample}_trimmed_R2.fq.gz \
        -j report_files/${sample}_fastp.json
        """
}

process runMultiQC {
        cpus {8}
        publishDir "${out_dir_qc}/", mode: 'copy', otherwise: false
        input:
                file('*') from json_report.collect()
        output:
                file("multiqc_report.html")
        """
        multiqc .
        """
}

