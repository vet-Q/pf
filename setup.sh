#!/bin/bash
# Pf Nextflow 파이프라인 - 환경 점검 스크립트 (Linux)

set -e

echo ""
echo "===================================================="
echo " Pf Nextflow 환경 점검"
echo "===================================================="
echo ""

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "프로젝트 디렉토리: $SCRIPT_DIR"
cd "$SCRIPT_DIR"

OK=0
FAIL=0

check_tool() {
    local name="$1"
    local cmd="$2"
    if command -v "$cmd" &>/dev/null; then
        local ver
        ver=$("$cmd" --version 2>&1 | head -1)
        echo "  ✅ $name  →  $ver"
        ((OK++))
    else
        echo "  ❌ $name  →  설치되어 있지 않습니다"
        ((FAIL++))
    fi
}

echo "[1/4] 필수 도구 확인"
echo "----------------------------------------------------"

# Java
if command -v java &>/dev/null; then
    JAVA_VER=$(java -version 2>&1 | head -1)
    echo "  ✅ Java  →  $JAVA_VER"
    ((OK++))
else
    echo "  ❌ Java  →  설치되어 있지 않습니다 (sudo apt install default-jdk)"
    ((FAIL++))
fi

# Nextflow
if command -v nextflow &>/dev/null; then
    NF_VER=$(nextflow -version 2>&1 | grep -oP 'version \K[0-9.]+' || echo "확인 불가")
    echo "  ✅ Nextflow  →  $NF_VER"
    ((OK++))
elif [ -x "$SCRIPT_DIR/nextflow" ]; then
    NF_VER=$("$SCRIPT_DIR/nextflow" -version 2>&1 | grep -oP 'version \K[0-9.]+' || echo "확인 불가")
    echo "  ✅ Nextflow (로컬)  →  $NF_VER"
    ((OK++))
else
    echo "  ❌ Nextflow  →  설치되어 있지 않습니다"
    echo "     설치: curl -s https://get.nextflow.io | bash && sudo mv nextflow /usr/local/bin/"
    ((FAIL++))
fi

check_tool "samtools" "samtools"
check_tool "minimap2" "minimap2"
check_tool "bcftools" "bcftools"
check_tool "ivar" "ivar"

# Rscript
if command -v Rscript &>/dev/null; then
    R_VER=$(Rscript --version 2>&1 | head -1)
    echo "  ✅ Rscript  →  $R_VER"
    ((OK++))
else
    echo "  ❌ Rscript  →  설치되어 있지 않습니다"
    ((FAIL++))
fi

echo ""
echo "[2/4] 참조 파일 확인"
echo "----------------------------------------------------"
REF="$SCRIPT_DIR/resources/PlasmoDB-67_Pfalciparum3D7_Genome.fasta"
BED="$SCRIPT_DIR/resources/amplicons.bed"

if [ -f "$REF" ]; then
    echo "  ✅ 레퍼런스 게놈: $(basename $REF)"
    if [ -f "${REF}.fai" ]; then
        echo "  ✅ FASTA 인덱스 (.fai) 존재"
    else
        echo "  ⚠️  FASTA 인덱스 없음 — samtools faidx 로 생성합니다..."
        samtools faidx "$REF" && echo "     → 인덱스 생성 완료" || echo "     → 인덱스 생성 실패"
    fi
else
    echo "  ❌ 레퍼런스 게놈 없음: $REF"
    ((FAIL++))
fi

if [ -f "$BED" ]; then
    echo "  ✅ BED 파일: $(basename $BED)"
else
    echo "  ❌ BED 파일 없음: $BED"
    ((FAIL++))
fi

echo ""
echo "[3/4] 입력 데이터 확인"
echo "----------------------------------------------------"
FASTQ_COUNT=$(find "$SCRIPT_DIR/data" -name "*.fastq.gz" -o -name "*.fastq" 2>/dev/null | wc -l)
BARCODE_COUNT=$(find "$SCRIPT_DIR/data" -mindepth 1 -maxdepth 1 -type d -name "barcode*" 2>/dev/null | wc -l)
echo "  📂 바코드 디렉토리: ${BARCODE_COUNT}개"
echo "  📄 FASTQ 파일 수: ${FASTQ_COUNT}개"

if [ "$FASTQ_COUNT" -eq 0 ]; then
    echo "  ⚠️  FASTQ 파일이 없습니다. data/barcode*/ 에 데이터를 넣어주세요."
fi

echo ""
echo "[4/4] bin 스크립트 실행 권한 확인"
echo "----------------------------------------------------"
NON_EXEC=0
for f in "$SCRIPT_DIR/bin/"*.{sh,R}; do
    [ -f "$f" ] || continue
    if [ ! -x "$f" ]; then
        echo "  ⚠️  실행 권한 없음: $(basename $f)"
        ((NON_EXEC++))
    fi
done

if [ "$NON_EXEC" -gt 0 ]; then
    echo "  → chmod +x bin/*.sh bin/*.R 로 권한을 부여합니다..."
    chmod +x "$SCRIPT_DIR/bin/"*.sh "$SCRIPT_DIR/bin/"*.R 2>/dev/null || true
    echo "  → 완료"
else
    echo "  ✅ 모든 스크립트에 실행 권한이 있습니다"
fi

echo ""
echo "===================================================="
if [ "$FAIL" -eq 0 ]; then
    echo " ✅ 환경 점검 완료  (${OK}개 도구 확인, 문제 없음)"
    echo ""
    echo " 파이프라인 실행:"
    echo "   ./run.sh"
    echo " 또는"
    echo "   nextflow run main.nf -profile standard"
else
    echo " ⚠️  환경 점검 완료  (${OK}개 통과 / ${FAIL}개 문제)"
    echo " 위의 ❌ 항목을 먼저 해결해 주세요."
fi
echo "===================================================="
echo ""
