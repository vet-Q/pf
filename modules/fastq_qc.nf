/*
 * FastQ Quality Control Module
 * Purpose: Assess quality of raw sequencing reads using FastQC
 * Input: Raw FASTQ files
 * Output: FastQC reports (HTML and zip)
 */

process FASTQ_QC {
    tag "${sample_id}"
    label "fastqc"
    publishDir "${params.outdir}/qc/fastqc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(fastq_files)
    
    output:
    path "*.html", emit: html_reports
    path "*_fastqc.zip", emit: zip_reports
    tuple val(sample_id), path(fastq_files), emit: fastq_pass  // Pass through for next step
    
    script:
    """
    # Run FastQC on all input files
    fastqc --outdir . --threads ${task.cpus} ${fastq_files.join(' ')}
    """
}

workflow FASTQ_QC_WF {
    take:
    fastq_ch  // Channel with tuple(sample_id, [fastq_files])
    
    main:
    FASTQ_QC(fastq_ch)
    
    emit:
    html_reports = FASTQ_QC.out.html_reports
    zip_reports = FASTQ_QC.out.zip_reports
    fastq_pass = FASTQ_QC.out.fastq_pass
}
