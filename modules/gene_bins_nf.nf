process GENE_BIN_DEPTH {

  tag { barcode }
  publishDir "${outdir}/gene_bins", mode: 'copy', overwrite: true

  input:
    tuple val(barcode), path(bam)
    path bed
    val outdir

  output:
    tuple val(barcode), path("${barcode}"), val(outdir)

  script:
  """
  set -euo pipefail

  echo "[DEBUG] PWD = \$(pwd)"
  echo "[DEBUG] bed path = ${bed}"
  ls -l ${bed}

  # Save absolute paths before changing directories
  bed_path=\$(readlink -f ${bed})
  bam_path=\$(readlink -f ${bam})

  mkdir -p ${barcode}
  cd ${barcode}

  while read -r chr start end gene; do
    [ -z "\$chr" ] && continue
    [ -z "\$gene" ] && continue

    region="\${chr}:\${start}-\${end}"
    echo "▶ ${barcode}  \$gene  \$region"

    mkdir -p "\$gene"
    cd "\$gene"

    samtools view -h "\$bam_path" "\$region" | \
    awk '
    BEGIN { OFS="\\t" }
    /^@/ {
      print > "bin_0-1k.sam"
      print > "bin_1-2k.sam"
      print > "bin_2-3k.sam"
      print > "bin_3-4k.sam"
      print > "bin_4kplus.sam"
      next
    }
    {
      len = length(\$10)
      if      (len < 1000)   print >> "bin_0-1k.sam";
      else if (len < 2000)   print >> "bin_1-2k.sam";
      else if (len < 3000)   print >> "bin_2-3k.sam";
      else if (len < 4000)   print >> "bin_3-4k.sam";
      else                   print >> "bin_4kplus.sam";
    }'

    for sam in bin_*.sam; do
      bam_out="\${sam%.sam}.bam"
      samtools view -bS "\$sam" | samtools sort -o "\$bam_out"
      samtools index "\$bam_out"
    done

    for b in bin_*.bam; do
      d="\${b%.bam}_depth.txt"
      samtools depth -aa -r "\$region" "\$b" > "\$d"
    done

    samtools view "\$bam_path" "\$region" | awk '{print length(\$10)}' > readlen.txt

    cd ..
  done < "\$bed_path"
  """
}

