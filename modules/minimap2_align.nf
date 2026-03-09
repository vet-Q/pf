/*
 * Minimap2 Read Alignment Module
 * Purpose: Align fastq reads to reference genome using minimap2
 * Input: Raw FASTQ files and reference genome
 * Output: Sorted and indexed BAM files
 */

process MINIMAP2_ALIGN {
    tag "${sample_id}"
    label "minimap2"
    publishDir "${params.outdir}/alignment/minimap2", mode: 'copy', pattern: "*.bam*"
    
    input:
    tuple val(sample_id), path(fastq_files)
    path reference
    val preset
    
    output:
    tuple val(sample_id), path("${sample_id}.Pf3D7.sorted.bam"), emit: bam
    
    script:
    """
    # Concatenate all FASTQ files for this sample, then align
    cat ${fastq_files.join(' ')} > combined_reads.fastq.gz

    minimap2 -a -x ${preset} \
        -t ${task.cpus} \
        "${reference}" \
        combined_reads.fastq.gz \
        > "${sample_id}.minimap2.sam"
    
    # Convert SAM to BAM and sort
    samtools view -@ ${task.cpus} -b "${sample_id}.minimap2.sam" | \
    samtools sort -@ ${task.cpus} -o "${sample_id}.Pf3D7.sorted.bam"
    
    # Index BAM file
    samtools index "${sample_id}.Pf3D7.sorted.bam"

    # Clean up intermediates
    rm -f combined_reads.fastq.gz "${sample_id}.minimap2.sam"
    """
}

process MINIMAP2_ALIGN_PRESET {
    tag "${sample_id}"
    label "minimap2"
    publishDir "${params.outdir}/alignment/minimap2", mode: 'copy', pattern: "*.bam*"
    
    input:
    tuple val(sample_id), path(fastq_files)
    path reference
    val preset  // "map-pb" (PacBio), "map-ont" (Oxford Nanopore), "sr" (short reads), "asm5" (asm)
    
    output:
    tuple val(sample_id), path("${sample_id}.Pf3D7.sorted.bam"), path("${sample_id}.Pf3D7.sorted.bam.bai"), emit: bam
    path "*.minimap2.sam", emit: sam
    
    script:
    """
    # Determine input format based on file extension
    minimap2 -a -x ${preset} \
        -t ${task.cpus} \
        "${reference}" \
        ${fastq_files.join(' ')} \
        > "${sample_id}.minimap2.sam"
    
    # Convert SAM to BAM and sort
    samtools view -@ ${task.cpus} -b "${sample_id}.minimap2.sam" | \\
    samtools sort -@ ${task.cpus} -o "${sample_id}.Pf3D7.sorted.bam"
    
    # Index BAM file
    samtools index "${sample_id}.Pf3D7.sorted.bam"
    """
}

workflow MINIMAP2_ALIGN_WF {
    take:
    fastq_ch    // Channel with tuple(sample_id, [fastq_files])
    reference   // Path to reference genome
    preset      // Alignment preset (optional, default: "sr")
    
    main:
    if (preset != null) {
        MINIMAP2_ALIGN_PRESET(fastq_ch, reference, preset)
        bam_ch = MINIMAP2_ALIGN_PRESET.out.bam
    } else {
        MINIMAP2_ALIGN(fastq_ch, reference)
        bam_ch = MINIMAP2_ALIGN.out.bam
    }
    
    emit:
    bam = bam_ch
}
