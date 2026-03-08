#!/bin/bash
# pfNextflow 초기 설정 스크립트 (Git Bash / WSL2)

set -e

echo ""
echo "===================================================="
echo " pfNextflow 초기 설정 스크립트"
echo "===================================================="
echo ""

# 현재 디렉토리 확인
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "[1/3] 프로젝트 디렉토리: $SCRIPT_DIR"
cd "$SCRIPT_DIR"

# Java 설정
echo ""
echo "[2/3] Java 설정..."
export JAVA_HOME="/c/Program Files/Java/jdk-25.0.2"
export PATH="$JAVA_HOME/bin:$PATH"

if "$JAVA_HOME/bin/java" -version 2>&1 | grep -q "version"; then
    echo "  ✓ Java 설치됨:"
    "$JAVA_HOME/bin/java" -version
else
    echo "  ✗ Java를 찾을 수 없습니다"
    echo "  설치 가이드: NEXTFLOW_INSTALL.md 참고"
    exit 1
fi

# Nextflow 확인 및 설정
echo ""
echo "[3/3] Nextflow 확인..."

if [ ! -f "nextflow" ]; then
    echo "  Nextflow 다운로드 중..."
    curl -s https://get.nextflow.io | bash
    chmod +x nextflow
fi

echo "  ✓ Nextflow 버전:"
./nextflow -version

echo ""
echo "===================================================="
echo " ✓ 설정 완료!"
echo "===================================================="
echo ""
echo "다음 명령으로 파이프라인 실행:"
echo "  nextflow run main.nf"
echo ""
echo "또는 다음 명령을 .bash_profile에 추가한 후 매번 실행:"
echo "  export JAVA_HOME='/c/Program Files/Java/jdk-25.0.2'"
echo "  export PATH=\"\$JAVA_HOME/bin:\$PATH\""
echo ""
