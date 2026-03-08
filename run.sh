#!/bin/bash
# WSL에서 Nextflow 한 번에 설치 및 실행

cd "/c/Users/kwono/OneDrive/문서/python project/pfNextflow"

# Java 설정
export JAVA_HOME="/c/Program Files/Java/jdk-25.0.2"
export PATH="$JAVA_HOME/bin:$PATH"

echo "=========================================="
echo "  pfNextflow - WSL 간단 실행"
echo "=========================================="
echo ""

# Step 1: Java 확인
echo "[1/3] Java 확인..."
"$JAVA_HOME/bin/java" -version
echo ""

# Step 2: Nextflow 설치
echo "[2/3] Nextflow 설치..."
if [ ! -f "nextflow" ]; then
    curl -s https://get.nextflow.io | bash
    chmod +x nextflow
fi
./nextflow -version
echo ""

# Step 3: 파이프라인 실행
echo "[3/3] 파이프라인 실행 중..."
echo ""
./nextflow run main.nf

