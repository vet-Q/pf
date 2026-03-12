/*
 * Minimap2 Read Alignment Module
 * Input : tuple(sample_id, [fastq_files]), reference FASTA, preset string
 * Output: tuple(sample_id, sorted.bam)
 */

process MINIMAP2_ALIGN {
    tag { "${sample_id}" }
    publishDir "${params.outdir}/alignment/minimap2", mode: 'copy', pattern: "*.bam*"

    input:
    tuple val(sample_id), path(fastq_files)
    path reference
    val  preset

    output:
    tuple val(sample_id), path("${sample_id}.Pf3D7.sorted.bam"), emit: bam

    script:
    """
    set -euo pipefail

    # 여러 FASTQ 파일을 하나로 합쳐서 정렬
    cat ${fastq_files.join(' ')} > combined_reads.fastq.gz

    minimap2 -a -x ${preset} \
        -t ${task.cpus} \
        "${reference}" \
        combined_reads.fastq.gz \
        > "${sample_id}.minimap2.sam"

    samtools view -@ ${task.cpus} -b "${sample_id}.minimap2.sam" | \
    samtools sort -@ ${task.cpus} -o "${sample_id}.Pf3D7.sorted.bam"

    samtools index "${sample_id}.Pf3D7.sorted.bam"

    rm -f combined_reads.fastq.gz "${sample_id}.minimap2.sam"
    """
}
