@echo off
REM Docker 이미지 사전 다운로드 스크립트
REM 이 스크립트는 선택사항입니다. 실행 시 약 5-10분 소요

setlocal enabledelayedexpansion

echo.
echo ================================================================================
echo  Docker 이미지 사전 다운로드
echo ================================================================================
echo.
echo 이 스크립트는 파이프라인에 필요한 모든 Docker 이미지를 미리 다운로드합니다.
echo 인터넷이 잘 연결되어 있고, 디스크 공간이 충분한지 확인하세요. (약 5GB)
echo.
pause

REM Docker 실행 확인
echo [1] Docker 확인...
docker --version 2>&1
if %errorlevel% neq 0 (
    echo ERROR: Docker를 찾을 수 없습니다!
    echo Docker Desktop을 설치하고 실행 중인지 확인하세요.
    pause
    exit /b 1
)

echo.
echo Docker 이미지 다운로드를 시작합니다...
echo.

REM 필요한 이미지 목록
setlocal
set images=^
    biocontainers/samtools:latest^
    biocontainers/fastqc:v0.12.1-1-deb_cv1^
    biocontainers/minimap2:v2.24-1-deb_cv1^
    andersenlabapps/ivar:latest^
    biocontainers/bcftools:v1.17-1-deb_cv1^
    biocontainers/bedtools:v2.31.0-1-deb_cv1^
    rocker/r-base:4.3.1

for %%i in (%images%) do (
    echo.
    echo [Downloading] %%i
    docker pull %%i
    if !errorlevel! neq 0 (
        echo WARNING: %%i 다운로드 실패 (네트워크 문제?)
    )
)

echo.
echo ================================================================================
echo  ✓ 다운로드 완료!
echo ================================================================================
echo.
echo 다운로드된 이미지 목록:
docker images | findstr /R "biocontainers andersenlabapps rocker"
echo.
echo 이제 run_pipeline.bat을 실행할 수 있습니다.
echo.
pause
