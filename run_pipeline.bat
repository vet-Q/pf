@echo off
REM pfNextflow 파이프라인 실행 스크립트

setlocal enabledelayedexpansion

echo.
echo ================================================================================
echo  pfNextflow 파이프라인 실행
echo ================================================================================
echo.

REM 현재 디렉토리
cd /d "%~dp0"

REM Java 재설정 (혹시 모르니)
echo [1] Java 설정...
set "JAVA_HOME=C:\Program Files\Java\jdk-25.0.2"
set "Path=!JAVA_HOME!\bin;!Path!"

java -version 2>&1
if %errorlevel% neq 0 (
    echo.
    echo ERROR: Java를 찾을 수 없습니다!
    echo 다시 시도:
    echo   1. 현재 모든 명령 프롬프트/PowerShell 닫기
    echo   2. 시스템 재부팅
    echo   3. 다시 이 파일 실행
    echo.
    pause
    exit /b 1
)

echo.
echo [2] Nextflow 확인...
if not exist "nextflow.cmd" (
    echo ERROR: nextflow.cmd를 찾을 수 없습니다!
    echo 먼저 install_nextflow.bat을 실행하세요.
    echo.
    pause
    exit /b 1
)

echo.
echo [3] configs/params.yaml 확인...
if not exist "configs\params.yaml" (
    echo ERROR: configs/params.yaml를 찾을 수 없습니다!
    echo.
    pause
    exit /b 1
)

echo.
echo ================================================================================
echo  파이프라인 시작
echo ================================================================================
echo.
echo 설정 파일: configs/params.yaml
echo 출력 디렉토리: (params.yaml의 paths.outdir 참고)
echo.

REM Nextflow 실행
call nextflow.cmd run main.nf -profile docker

if %errorlevel% equ 0 (
    echo.
    echo ================================================================================
    echo  ✓ 파이프라인 완료!
    echo ================================================================================
) else (
    echo.
    echo ================================================================================
    echo  ✗ 파이프라인 오류 발생
    echo ================================================================================
    echo.
    echo 문제 해결:
    echo   1. Docker Desktop이 실행 중인지 확인
    echo   2. configs/params.yaml 설정 확인
    echo   3. 로그 파일 확인: .nextflow.log
    echo.
)

pause
