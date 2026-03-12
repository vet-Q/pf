#!/bin/bash
# Plasmodium falciparum Nextflow Pipeline - Linux Runner

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Auto-detect Java home (Windows or Linux)
if [ -d "/usr/lib/jvm/java-21-openjdk-amd64" ]; then
    # Linux
    export JAVA_HOME="/usr/lib/jvm/java-21-openjdk-amd64"
elif [ -d "/c/java/jdk-25.0.2" ]; then
    # Windows (copied Java)
    export JAVA_HOME="/c/java/jdk-25.0.2"
elif [ -d "/c/Program Files/Java/jdk-25.0.2" ]; then
    # Windows (original path)
    export JAVA_HOME="/c/Program Files/Java/jdk-25.0.2"
fi
export PATH="$JAVA_HOME/bin:$PATH"

echo "=========================================="
echo "  pfNextflow - Linux Pipeline Runner"
echo "=========================================="
echo ""

# Step 1: Check Java
echo "[1/4] Checking Java..."
if ! command -v java &>/dev/null; then
    echo "ERROR: Java not found. Install Java 11+ (e.g. sudo apt install default-jdk)"
    exit 1
fi
java -version 2>&1 | head -1
echo ""

# Step 2: Check Nextflow
echo "[2/4] Checking Nextflow..."
if [ -f "./nextflow" ]; then
    NF="./nextflow"
elif command -v nextflow &>/dev/null; then
    NF="nextflow"
else
    echo "Nextflow not found. Installing..."
    curl -s https://get.nextflow.io | bash
    chmod +x nextflow
    NF="./nextflow"
fi
$NF -version
echo ""

# Step 3: Check required tools
echo "[3/4] Checking bioinformatics tools..."
for tool in samtools minimap2 bcftools ivar mafft Rscript; do
    if command -v "$tool" &>/dev/null; then
        echo "  ✓ $tool found: $(which $tool)"
    else
        echo "  ✗ $tool NOT found - please install it"
    fi
done
echo ""

# Step 4: Index reference genome if needed
REF="resources/PlasmoDB-67_Pfalciparum3D7_Genome.fasta"
if [ -f "$REF" ] && [ ! -f "${REF}.fai" ]; then
    echo "  Indexing reference genome..."
    samtools faidx "$REF"
fi
echo ""

# Step 5: Run pipeline (local profile, no Docker)
echo "[4/4] Running pipeline..."
echo ""
$NF run main.nf -profile standard

