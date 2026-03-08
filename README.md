# Plasmodium falciparum (Pf) Analysis Pipeline

A comprehensive Nextflow-based pipeline for analyzing Plasmodium falciparum genomic data. This pipeline processes raw FASTQ files or pre-aligned BAM files to identify variants, generate consensus sequences, analyze coverage depths, and produce publication-ready visualizations focused on drug-resistance genes.

## 🚀 빠른 시작 (5분)

**복잡함이 없이 시작하고 싶으신가요?**

👉 **[SIMPLE_START.md](SIMPLE_START.md)** ← 여기부터 시작하세요!

간단한 3단계만으로 끝납니다:
1. D:\input 폴더 만들기 + FASTQ 파일 넣기
2. `run.bat D:\input D:\output conda` 실행하기
3. 결과 확인 (D:\output)

**CMD/PowerShell 창에서 한줄이면 끝:**
```cmd
cd C:\Users\kwono\OneDrive\문서\python project\pfNextflow
run.bat D:\input D:\output conda
```

---

## Key Features

✅ **Complete from raw reads**: FASTQ QC → Alignment → Variant calling → Visualization
✅ **Docker-based (Offline capable)**: All tools containerized, runs without external dependencies
✅ **Flexible input**: Accept both raw FASTQ and pre-aligned BAM files
✅ **Drug-resistance profiling**: Focused analysis on key Pf resistance genes
✅ **Publication-ready plots**: Comprehensive gene-level visualizations

## Overview

This pipeline is designed for research groups analyzing Pf samples, particularly those interested in:
- **FastQ quality assessment** using FastQC
- **Read alignment** from fastq files using minimap2 (supports PacBio, Nanopore, short reads)
- **Variant calling** from aligned sequencing data with bcftools
- **Consensus sequence** generation using iVar
- **Coverage depth analysis** across gene regions
- **Gene-specific visualizations** including allele frequency and depth plots
- **Drug-resistance gene profiling** (CRT, DHFR, DHPS, K13, MDR1, MSP2, PMI, PMIII)

## Requirements

### System Requirements
- **Nextflow** >= 22.0
- **Docker** (recommended for offline operation) OR **Conda** (for local environment)
- **Disk space** >= 100GB (for typical multi-sample analyses)
- **RAM** >= 16GB (configurable)

### Option A: Docker (RECOMMENDED - Offline Capable)
All bioinformatics tools are containerized and managed by Nextflow automatically.
```bash
# Install Docker Desktop
# macOS/Windows: https://www.docker.com/products/docker-desktop
# Linux: sudo apt-get install docker.io (or your package manager equivalent)

# Verify installation
docker --version
```

### Option B: Conda (Alternative)
For local installation without Docker:
```bash
# Install Miniconda
# https://docs.conda.io/projects/conda/en/latest/user-guide/install/

# Verify installation
conda --version
```

### Bioinformatics Tools (Managed Automatically)
The pipeline relies on the following tools (all automatically managed):
- `fastqc` - FastQ quality control
- `minimap2` - Long-read and short-read alignment
- `samtools` - BAM file processing
- `bcftools` - Variant calling and VCF manipulation
- `iVar` - Consensus sequence generation
- `bedtools` - Genomic region operations
- `R` (>= 4.0) - Visualization and statistical analysis

## Installation

### 1. Clone the Repository
```bash
git clone https://github.com/yourusername/pfNextflow.git
cd pfNextflow
```

### 2. Install Nextflow
```bash
# Using curl
curl -s https://get.nextflow.io | bash
chmod +x nextflow

# Or using Homebrew (macOS/Linux)
brew install nextflow
```

### 3. Verify Installation
```bash
./nextflow version
```

## Quick Start

### Basic Usage

```bash
# Run with default parameters (requires input BAM files)
nextflow run main.nf \
  --bam_root /path/to/bam/directory \
  --outdir /path/to/output \
  -profile docker

# Run with custom gene list
nextflow run main.nf \
  --bam_root /path/to/bam/directory \
  --outdir /path/to/output \
  --genes "CRT DHFR K13 MDR1" \
  -profile docker
```

### Example: Running All Analyses
```bash
nextflow run main.nf \
  --bam_root /data/bams \
  --outdir /data/results \
  --do_depth true \
  --do_ivar true \
  --do_variants true \
  --do_multipanel true \
  --do_gene_1x8 true \
  --do_afdp true \
  -profile docker
```

### Example: From Raw FASTQ Files
```bash
nextflow run main.nf \
  --fastq_root /data/fastq \
  --outdir /data/results \
  --do_fastqc true \
  --do_alignment true \
  --minimap2_preset "sr" \
  --do_variants true \
  --do_afdp true \
  -profile docker
```

### Example: Nanopore Long Reads
```bash
nextflow run main.nf \
  --fastq_root /data/fastq \
  --outdir /data/results \
  --do_fastqc true \
  --do_alignment true \
  --minimap2_preset "map-ont" \
  --do_variants true \
  -profile docker
```

## Docker Offline Usage

### Pre-download Docker Images (Offline Capability)
```bash
# Pull all required images to local disk
docker pull biocontainers/fastqc:v0.12.1-1-deb_cv1
docker pull biocontainers/minimap2:v2.24-1-deb_cv1
docker pull biocontainers/bcftools:v1.17-1-deb_cv1
docker pull biocontainers/bedtools:v2.31.0-1-deb_cv1
docker pull andersenlabapps/ivar:latest
docker pull rocker/r-base:4.3.1
docker pull biocontainers/samtools:latest

# Save to tar files for transfer
mkdir -p docker-images
docker save biocontainers/fastqc:v0.12.1-1-deb_cv1 | gzip > docker-images/fastqc.tar.gz
docker save biocontainers/minimap2:v2.24-1-deb_cv1 | gzip > docker-images/minimap2.tar.gz
# ... repeat for all images
```

### Load Pre-downloaded Images on Offline Machine
```bash
# Load all images
for file in docker-images/*.tar.gz; do
    docker load < "$file"
done

# Verify images are loaded
docker images
```

### Run Pipeline Offline
```bash
nextflow run main.nf \
  --fastq_root /data/fastq \
  --outdir /data/results \
  --do_fastqc true \
  --do_alignment true \
  --do_variants true \
  -profile docker \
  -offline
```

## Pipeline Parameters

### Input Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `bam_root` | `.` | Root directory for BAM files (when not aligning) |
| `bam_glob` | `merge/nomadic/minknow/barcodes/*/bams/*.Pf3D7.final.sorted.bam` | Glob pattern to find BAM files |
| `fastq_root` | `.` | Root directory for FASTQ files (for alignment) |
| `fastq_glob` | `**/*.{fastq,fq,fastq.gz,fq.gz}` | Glob pattern to find FASTQ files |
| `bed` | `resources/amplicons.bed` | BED file with amplicon coordinates |
| `ref` | `resources/PlasmoDB-67_Pfalciparum3D7_Genome.fasta` | Reference genome FASTA |
| `outdir` | `../data/results/expX` | Output directory path |
| `genes` | `CRT DHFR DHPS K13 MDR1 MSP2 PMI PMIII` | Space-separated gene names for analysis |

### FastQ Quality Control & Alignment

| Parameter | Default | Description |
|-----------|---------|-------------|
| `do_fastqc` | `false` | Run FastQC quality assessment on raw reads |
| `do_alignment` | `false` | Perform read alignment using minimap2 |
| `minimap2_preset` | `sr` | Alignment preset: "sr" (short reads), "map-pb" (PacBio), "map-ont" (Nanopore), "asm5" (assembly) |
| `minimap2_threads` | `4` | Number of threads for minimap2 |

### iVar Consensus Generation

| Parameter | Default | Description |
|-----------|---------|-------------|
| `ivar_threshold` | `0.6` | Minimum allele frequency threshold (0-1) |
| `ivar_min_cons_cov` | `10` | Minimum depth for consensus calling |
| `ivar_min_var_cov` | `1` | Minimum depth for variant calling |

### Variant Calling (bcftools)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `variant_method` | `bcftools` | Variant calling method |
| `bcftools_mode` | `mv` | bcftools mode ("m"=multiallelic, "v"=vcf_outputs) |
| `snp_only` | `true` | Only report SNPs (exclude indels) |
| `bcftools_threads` | `4` | Number of threads for bcftools |
| `max_depth` | `10000` | Maximum depth threshold (filters high-depth regions) |

### AF/DP Figure Generation

| Parameter | Default | Description |
|-----------|---------|-------------|
| `afdp_min_depth` | `15` | Minimum depth for plotting |
| `afdp_af_threshold` | `0.50` | Allele frequency threshold (0-1) |
| `afdp_shade_xmin` | `0.48` | Shading region minimum (x-axis) |
| `afdp_shade_xmax` | `0.52` | Shading region maximum (x-axis) |

### Analysis Toggles (Enable/Disable)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `do_fastqc` | `false` | Perform quality assessment on raw reads |
| `do_alignment` | `false` | Align reads using minimap2 (requires FASTQ input) |
| `do_depth` | `false` | Run depth analysis |
| `do_ivar` | `false` | Generate consensus sequences |
| `do_gene_bins` | `false` | Calculate gene bin depth (required for plots) |
| `do_variants` | `false` | Call variants with bcftools |
| `do_multipanel` | `false` | Generate multipanel gene plots |
| `do_gene_1x8` | `false` | Generate 1x8 gene comparison plots |
| `do_afdp` | `true` | Generate AF/DP figures |

## Output Structure

```
results/
├── qc/
│   └── fastqc/
│       ├── sample1_fastqc.html
│       ├── sample1_fastqc.zip
│       └── ...
├── alignment/
│   └── minimap2/
│       ├── sample1.Pf3D7.sorted.bam
│       ├── sample1.Pf3D7.sorted.bam.bai
│       ├── sample1.minimap2.sam
│       └── ...
├── depth/
│   ├── sample1/
│   │   └── depth_plots.pdf
│   └── sample2/
│       └── depth_plots.pdf
├── consensus/
│   ├── sample1.consensus.fasta
│   └── sample2.consensus.fasta
├── calling/
│   └── bcftools/
│       ├── sample1.vcf.gz
│       ├── sample1.vcf.gz.tbi
│       └── ...
├── gene_bins/
│   └── gene_bin_depth_data.csv
├── gene_plots/
│   ├── multipanel_genes.pdf
│   └── gene_1x8_comparison.pdf
├── afdp_figures/
│   └── allele_freq_depth.pdf
└── pipeline_info/
    ├── execution_timeline_*.html
    ├── execution_report_*.html
    ├── execution_trace_*.txt
    └── pipeline_dag_*.pdf
```

## Advanced Usage

### Execution Profiles

```bash
# Docker (RECOMMENDED - all tools containerized, offline capable)
nextflow run main.nf -profile docker

# Conda (local installation, requires all tools installed)
nextflow run main.nf -profile conda

# SLURM HPC Cluster
nextflow run main.nf -profile slurm

# Singularity (alternative to Docker)
nextflow run main.nf -profile singularity

# Local development/testing
nextflow run main.nf -profile local
```

### Running with Conda Profiles
```bash
nextflow run main.nf \
  --fastq_root /data/fastq \
  --outdir /data/results \
  --do_fastqc true \
  --do_alignment true \
  -profile conda
```

### Resuming Failed Runs
```bash
nextflow run main.nf \
  --fastq_root /data/fastq \
  --outdir /data/results \
  -resume
```

### Custom Reference Files
```bash
nextflow run main.nf \
  --fastq_root /data/fastq \
  --ref /custom/path/reference.fasta \
  --bed /custom/path/amplicons.bed \
  --outdir /data/results \
  -profile docker
```

### Parallel Execution with Reporting
```bash
nextflow run main.nf \
  --fastq_root /data/fastq \
  --outdir /data/results \
  -profile docker \
  -with-timeline timeline.html \
  -with-report report.html \
  -with-trace trace.txt
```

### Running with Custom Resource Limits
```bash
nextflow run main.nf \
  --fastq_root /data/fastq \
  --outdir /data/results \
  --max_cpus 16 \
  --max_memory 64.GB \
  -profile docker
```

## Project Structure

```
pfNextflow/
├── main.nf                          # Main workflow
├── nextflow.config                  # Nextflow configuration (Docker, Conda, profiles)
├── modules/                         # Analysis modules
│   ├── fastq_qc.nf                  # FastQC quality assessment
│   ├── minimap2_align.nf            # Read alignment (minimap2)
│   ├── pf_depth.nf                  # Depth analysis
│   ├── ivar_consensus.nf            # Consensus generation
│   ├── call_variants.nf             # Variant calling
│   ├── gene_bins_nf.nf              # Gene bin depth
│   ├── gene_bins_plot_nf.nf         # Multipanel plots
│   ├── gene_1to8_plot.nf            # 1x8 comparison
│   └── af_dp_figures.nf             # AF/DP visualizations
├── environments/                    # Conda environment specs
│   ├── default.yaml                 # Default environment (all tools)
│   ├── fastqc.yaml                  # FastQC environment
│   ├── minimap2.yaml                # Minimap2 + samtools
│   ├── ivar.yaml                    # iVar environment
│   ├── bcftools.yaml                # bcftools environment
│   └── r-analysis.yaml              # R visualization environment
├── bin/                             # Helper scripts
│   ├── call_bcftools.sh             # bcftools wrapper
│   ├── call_ivar.sh                 # iVar wrapper
│   ├── make_consensus_*.sh          # Consensus generation
│   ├── make_af_dp_figures.R         # R visualization
│   ├── plot_depth_*.R               # Depth plotting
│   └── plot_gene_*.R                # Gene-level plots
├── resources/                       # Reference files
│   ├── amplicons.bed                # Gene amplicon regions
│   ├── PlasmoDB-67_Pfalciparum3D7_Genome.fasta  # Reference genome
│   ├── metadata_country.csv         # Sample metadata
│   └── ...
├── .gitignore                       # Git ignore patterns
├── README.md                        # This file
└── test.sh                          # Testing script
```

## Data Preparation

### Input Requirements
1. **BAM files** must be:
   - Aligned to Pf3D7 reference genome
   - Sorted by coordinate
   - Indexed (.bai files present)
   - Named with pattern: `*.Pf3D7.final.sorted.bam`

2. **Reference files** (included):
   - BED file with amplicon coordinates
   - Pf3D7 reference genome FASTA

### Creating Index Files
```bash
# Index BAM files if needed
samtools index sample.sorted.bam

# Index reference FASTA if needed
samtools faidx reference.fasta
```

## Troubleshooting

### Common Issues

#### Issue: "No FASTQ files found"
```bash
# Check that FASTQ files exist and match the glob pattern
ls -la /your/fastq/root/**/*.fastq.gz

# Update parameters for custom structure
nextflow run main.nf --fastq_root /actual/path --fastq_glob "*.fastq.gz"
```

#### Issue: Docker daemon not running
```bash
# Start Docker Desktop or Docker daemon
# macOS: Open Docker app from Applications
# Linux: sudo systemctl start docker
# Windows: Start Docker Desktop

# Verify Docker is running
docker ps
```

#### Issue: Docker image download fails (Network issue)
```bash
# Use offline mode with pre-downloaded images
nextflow run main.nf -offline -profile docker

# Or load images directly
docker load < docker-images/fastqc.tar.gz
```

#### Issue: Permission denied in work directory
```bash
# Check work directory permissions
chmod -R u+w work/
nextflow clean -f
```

#### Issue: Out of memory during alignment
```bash
# Run with reduced resources
nextflow run main.nf \
  --fastq_root /data/fastq \
  --minimap2_threads 2 \
  --max_memory 8.GB \
  -profile docker
```

#### Issue: "No BAM files found" (for BAM input)
```bash
# Check that BAM files exist and match pattern
ls -la /your/bam/root/merge/nomadic/minknow/barcodes/*/bams/*.Pf3D7.final.sorted.bam

# Verify indexes exist
ls -la *.bai

# Update parameters if different structure
nextflow run main.nf --bam_root /actual/path --bam_glob "**/*.bam"
```

#### Issue: Missing reference files
```bash
# Verify reference files exist
ls -la resources/
ls -la resources/PlasmoDB-67_Pfalciparum3D7_Genome.fasta
samtools faidx resources/PlasmoDB-67_Pfalciparum3D7_Genome.fasta
```

#### Issue: Conda environment creation fails
```bash
# Update conda
conda update -n base conda

# Clean cache
conda clean --all

# Recreate environment
nextflow run main.nf -profile conda --clean
```

### Docker-Specific Troubleshooting

#### Check available images
```bash
docker images | grep -E "fastqc|minimap2|bcftools"
```

#### Remove stale containers/images
```bash
# Remove stopped containers
docker container prune

# Remove dangling images
docker image prune
```

#### View Docker execution logs
```bash
# Enable Nextflow debug mode
nextflow run main.nf -debug -profile docker
```

### Offline Mode Checklist

- [ ] All Docker images downloaded and saved locally
- [ ] Docker service running on target machine
- [ ] Images loaded into Docker daemon (`docker load`)
- [ ] Nextflow cache cleared (`nextflow clean -f`)
- [ ] Running with `-offline` flag
- [ ] All reference data available locally

## Performance Tips

1. **Enable caching**: Nextflow automatically caches intermediate results
2. **Use parallel processing**: Adjust `bcftools_threads` based on CPU availability
3. **Monitor progress**: Use Nextflow reports with `-with-timeline` and `-with-report`
4. **Batch samples**: Process multiple samples in single runs efficiently

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch (`git checkout -b feature/your-feature`)
3. Commit changes (`git commit -am 'Add your feature'`)
4. Push to branch (`git push origin feature/your-feature`)
5. Open a Pull Request

## Citation

If you use this pipeline in your research, please cite:

```
[Author et al. (Year). Pipeline name. Journal/Repository]
```

## License

This project is licensed under the MIT License - see LICENSE file for details.

## Contact & Support

- **Issues**: Please report bugs and feature requests on GitHub Issues
- **Email**: [your-email@institution.edu]
- **Documentation**: Additional help at `docs/` (if available)

## Acknowledgments

- Built with [Nextflow](https://www.nextflow.io/) - Data-driven computational pipelines for complex genome analysis
- Reference genome from [PlasmoDB](https://plasmodb.org/)
- Variant analysis using [bcftools](http://samtools.github.io/bcftools/)
- Consensus generation using [iVar](https://andersen-lab.github.io/ivar/html/)

## Future Enhancements

- [ ] Support for additional variant callers (GATK, Freebayes)
- [ ] Integration with population genetics analysis
- [ ] Automated report generation
- [ ] Web-based visualization portal
- [ ] Expanded drug-resistance gene profiles
- [ ] Integration with ClonalFrameML for phylogenetic analysis

## Configuration Files

### nextflow.config
```
location: ./nextflow.config
Contains:
- Docker image specifications for each process
- Conda environment definitions
- HPC/SLURM settings
- Resource limits (CPU, memory, time)
- Report generation settings
```

### Environment YAML Files
Located in `environments/` directory:
- `fastqc.yaml` - FastQC dependencies
- `minimap2.yaml` - Minimap2 and samtools
- `ivar.yaml` - iVar, samtools, bcftools
- `bcftools.yaml` - Bcftools, samtools, bedtools
- `r-analysis.yaml` - R and visualization packages
- `default.yaml` - All tools combined

## Docker Container References

The pipeline automatically pulls these images:

| Tool | Image | Version |
|------|-------|---------|
| FastQC | `biocontainers/fastqc` | v0.12.1 |
| minimap2 | `biocontainers/minimap2` | v2.24 |
| samtools | `biocontainers/samtools` | latest |
| bcftools | `biocontainers/bcftools` | v1.17 |
| bedtools | `biocontainers/bedtools` | v2.31.0 |
| iVar | `andersenlabapps/ivar` | latest |
| R | `rocker/r-base` | 4.3.1 |

### Pre-cache Docker Images for Offline Use
```bash
# Create directory for Docker images
mkdir -p docker-images

# Pull all images
docker pull biocontainers/fastqc:v0.12.1-1-deb_cv1
docker pull biocontainers/minimap2:v2.24-1-deb_cv1
docker pull biocontainers/samtools:latest
docker pull biocontainers/bcftools:v1.17-1-deb_cv1
docker pull biocontainers/bedtools:v2.31.0-1-deb_cv1
docker pull andersenlabapps/ivar:latest
docker pull rocker/r-base:4.3.1

# Export to tar.gz files
docker save biocontainers/fastqc:v0.12.1-1-deb_cv1 | gzip > docker-images/fastqc.tar.gz
docker save biocontainers/minimap2:v2.24-1-deb_cv1 | gzip > docker-images/minimap2.tar.gz
docker save biocontainers/samtools:latest | gzip > docker-images/samtools.tar.gz
docker save biocontainers/bcftools:v1.17-1-deb_cv1 | gzip > docker-images/bcftools.tar.gz
docker save biocontainers/bedtools:v2.31.0-1-deb_cv1 | gzip > docker-images/bedtools.tar.gz
docker save andersenlabapps/ivar:latest | gzip > docker-images/ivar.tar.gz
docker save rocker/r-base:4.3.1 | gzip > docker-images/r-base.tar.gz

# Transfer docker-images/ folder to offline machine
```

### Load Images on Offline Machine
```bash
# Load all images
for file in docker-images/*.tar.gz; do
    gunzip -c "$file" | docker load
done

# Verify images
docker images
```

## Future Enhancements

- [ ] Support for additional variant callers (GATK, Freebayes)
- [ ] Integration with population genetics analysis
- [ ] Automated report generation with interactive dashboards
- [ ] Web-based visualization portal
- [ ] Expanded drug-resistance gene profiles
- [ ] Integration with ClonalFrameML for phylogenetic analysis
- [ ] Snakemake alternative workflow
- [ ] GPU acceleration for alignment

---

**Last Updated**: March 2026  
**Pipeline Version**: 1.0
