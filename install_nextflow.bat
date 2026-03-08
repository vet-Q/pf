@echo off
REM Nextflow 간단 설치 스크립트 (Windows Command Prompt / PowerShell)

setlocal enabledelayedexpansion

echo.
echo ================================================================================
echo  Nextflow 설치 및 첫 실행
echo ================================================================================
echo.

REM 현재 디렉토리
cd /d "%~dp0"
echo 현재 디렉토리: %cd%
echo.

REM Java 버전 확인
echo [Step 1] Java 확인...
java -version 2>&1
if %errorlevel% neq 0 (
    echo.
    echo ERROR: Java를 찾을 수 없습니다!
    echo 해결방법:
    echo   1. 새로운 PowerShell/CMD 창을 열기 (현재 창 종료)
    echo   2. 재부팅
    echo   3. Java 수동 설치: https://www.oracle.com/java/technologies/downloads/
    echo.
    pause
    exit /b 1
)

echo.
echo [Step 2] Docker 확인...
docker --version 2>&1 || (
    echo WARNING: Docker를 찾을 수 없습니다.
    echo Docker Desktop 설치: https://www.docker.com/products/docker-desktop
    echo.
)

echo.
echo [Step 3] Nextflow 설치 중...

REM 이미 nextflow가 있는지 확인
if exist "nextflow.cmd" (
    echo  ✓ Nextflow 이미 설치됨
) else (
    echo  Nextflow 다운로드 중... (약 1-2분)
    powershell -NoProfile -Command "$ProgressPreference='SilentlyContinue'; Invoke-WebRequest -Uri 'https://get.nextflow.io' -OutFile 'nextflow.cmd'; Write-Host ' ✓ Nextflow 다운로드 완료'"
)

echo.
echo [Step 4] Nextflow 버전 확인...
call nextflow.cmd -version 2>&1 || (
    echo ERROR: Nextflow 실행 실패
    pause
    exit /b 1
)

echo.
echo ================================================================================
echo  ✓ 설치 완료!
echo ================================================================================
echo.
echo 다음 명령으로 파이프라인 실행:
echo   nextflow run main.nf
echo.
echo 또는 run_pipeline.bat 파일을 클릭하세요
echo.
pause
