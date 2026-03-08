@echo off
REM pfNextflow 초기 설정 스크립트 (Windows)
REM 이 파일을 관리자 권한으로 실행하세요

setlocal enabledelayedexpansion

echo.
echo ====================================================
echo  pfNextflow 초기 설정
echo ====================================================
echo.

REM Java 확인
echo [1/3] Java 확인 중...
java -version >nul 2>&1
if %errorlevel% equ 0 (
    echo     ✓ Java 설치됨
    for /f "tokens=3" %%a in ('java -version 2^>^&1 ^| findstr "version"') do echo     버전: %%a
) else (
    echo     ✗ Java 없음 (이미 설치되어야 함)
    pause
    exit /b 1
)

REM Docker 확인
echo.
echo [2/3] Docker 확인 중...
docker --version >nul 2>&1
if %errorlevel% equ 0 (
    echo     ✓ Docker 설치됨
    docker --version
) else (
    echo     ✗ Docker 없음
    echo     Docker Desktop: https://www.docker.com/products/docker-desktop
    pause
    exit /b 1
)

REM Nextflow 확인 및 설치
echo.
echo [3/3] Nextflow 설치 중...

cd /d "%~dp0"

if exist "nextflow.cmd" (
    echo     ✓ Nextflow 이미 설치됨
    nextflow -version
) else (
    echo     Nextflow 다운로드 중...
    powershell -NoProfile -Command "$ProgressPreference='SilentlyContinue'; Invoke-WebRequest -Uri 'https://get.nextflow.io' -OutFile 'nextflow.cmd'"
    
    if exist "nextflow.cmd" (
        echo     ✓ Nextflow 설치 완료
        nextflow -version
    ) else (
        echo     ✗ Nextflow 다운로드 실패
        pause
        exit /b 1
    )
)

echo.
echo ====================================================
echo  ✓ 설정 완료!
echo ====================================================
echo.
echo 다음 명령으로 파이프라인을 실행할 수 있습니다:
echo   nextflow run main.nf
echo.
echo 더 많은 설정 옵션: configs\params.yaml 참고
echo.
pause
