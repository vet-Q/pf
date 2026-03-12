/*
 * FastQ Quality Control Module
 * Input : tuple(sample_id, [fastq_files])
 * Output: FastQC HTML/ZIP reports + pass-through fastq_pass tuple
 */

process FASTQ_QC {
    tag { "${sample_id}" }
    publishDir "${params.outdir}/qc/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq_files)

    output:
    path "*.html",         emit: html_reports
    path "*_fastqc.zip",   emit: zip_reports
    tuple val(sample_id), path(fastq_files), emit: fastq_pass

    script:
    """
    fastqc --outdir . --threads ${task.cpus} ${fastq_files.join(' ')}
    """
}
