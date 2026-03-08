# 🚀 Nextflow 설치 및 Docker 오프라인 실행 가이드

## 📋 시스템 요구사항

- **Java 11+** (필수)
- **Docker** (오프라인 모드로 실행)
- 또는 **WSL2** (Windows의 경우)

---

## 1️⃣ Step 1: Java 설치

### Windows에서 Java 설치

#### 방법 1️⃣: 자동 설치 (권장)

**Chocolatey를 사용하는 경우** (이미 설치되어 있다면):
```powershell
# PowerShell을 관리자 권한으로 열기
choco install openjdk11 -y
```

**또는 직접 다운로드:**
1. [Oracle JDK Download](https://www.oracle.com/kr/java/technologies/downloads/) 방문
2. **Java 17** 또는 **OpenJDK 11** 다운로드
3. 설치 후 **재부팅**

#### 설치 확인

```powershell
java -version
```

정상이면 다음과 같이 출력됩니다:
```
openjdk version "11.0.15" 2022-04-19
OpenJDK Runtime Environment
```

---

## 2️⃣ Step 2: Docker 설치

### Windows에서 Docker 설치

#### 방법 1️⃣: Docker Desktop (권장)

1. [Docker Desktop Download](https://www.docker.com/products/docker-desktop) 방문
2. **Windows** 버전 다운로드 (WSL2 필요)
3. 설치 실행
4. 재부팅

#### 방법 2️⃣: Chocolatey 사용

```powershell
# PowerShell 관리자 권한
choco install docker-desktop -y
```

#### 설치 확인

```powershell
docker --version
docker run hello-world
```

⚠️ **중요**: 첫 실행 시 Docker 데몬이 시작되어야 합니다. Docker Desktop 앱을 열어두세요.

---

## 3️⃣ Step 3: Nextflow 설치

### Windows에서 Nextflow 설치

#### 방법 1️⃣: Git Bash / WSL2에서 설치 (권장)

**WSL2 또는 Git Bash 열기:**

```bash
# 1️⃣ Nextflow 다운로드
curl -s https://get.nextflow.io | bash

# 2️⃣ nextflow 파일을 PATH에 추가
sudo mv nextflow /usr/local/bin/
chmod +x /usr/local/bin/nextflow

# 3️⃣ 설치 확인
nextflow -version
```

#### 방법 2️⃣: PowerShell에서 설치

```powershell
# 1️⃣ Nextflow 다운로드 URL에서 직접 다운로드
Invoke-WebRequest -Uri "https://github.com/nextflow-io/nextflow/releases/download/v24.04.4/nextflow-24.04.4-dist.tar.gz" -OutFile "nextflow.tar.gz"

# 2️⃣ 압축 해제 (7-Zip 또는 WinRAR 필요)
# 또는 WSL에서 해제하는 것이 더 간단

# 3️⃣ 환경 변수에 추가 (선택)
# "nextflow.bat"를 C:\Program Files\Nextflow\ 등에 배치
```

#### 방법 3️⃣: Chocolatey 사용

```powershell
choco install nextflow -y
```

#### 설치 확인

```bash
# WSL2 / Git Bash에서:
nextflow -version

# 정상 출력:
# nextflow version 24.04.4
```

---

## 4️⃣ Step 4: 오프라인 Docker 이미지 준비

파이프라인에 필요한 Docker 이미지들을 **미리 다운로드**해야 합니다.

### 필요한 이미지 목록

```bash
# 이미지 다운로드 (인터넷 연결 필요)
docker pull biocontainers/samtools:latest
docker pull biocontainers/fastqc:v0.12.1-1-deb_cv1
docker pull biocontainers/minimap2:v2.24-1-deb_cv1
docker pull andersenlabapps/ivar:latest
docker pull biocontainers/bcftools:v1.17-1-deb_cv1
docker pull biocontainers/bedtools:v2.31.0-1-deb_cv1
docker pull rocker/r-base:4.3.1
```

### 오프라인으로 이미지 저장/로드하기

**인터넷이 있는 컴퓨터에서:**

```bash
# 1️⃣ 이미지 저장
mkdir -p docker_images

docker save biocontainers/samtools:latest -o docker_images/samtools.tar
docker save biocontainers/fastqc:v0.12.1-1-deb_cv1 -o docker_images/fastqc.tar
docker save biocontainers/minimap2:v2.24-1-deb_cv1 -o docker_images/minimap2.tar
docker save andersenlabapps/ivar:latest -o docker_images/ivar.tar
docker save biocontainers/bcftools:v1.17-1-deb_cv1 -o docker_images/bcftools.tar
docker save biocontainers/bedtools:v2.31.0-1-deb_cv1 -o docker_images/bedtools.tar
docker save rocker/r-base:4.3.1 -o docker_images/r-base.tar

# 2️⃣ USB 또는 외부 저장소로 복사
# docker_images/ 폴더를 복사하세요
```

**오프라인 컴퓨터에서:**

```bash
# 1️⃣ USB에서 이미지 폴더 복사
# docker_images/ 폴더를 현재 디렉토리에 복사

# 2️⃣ 이미지 로드
docker load -i docker_images/samtools.tar
docker load -i docker_images/fastqc.tar
docker load -i docker_images/minimap2.tar
docker load -i docker_images/ivar.tar
docker load -i docker_images/bcftools.tar
docker load -i docker_images/bedtools.tar
docker load -i docker_images/r-base.tar

# 3️⃣ 로드된 이미지 확인
docker images
```

---

## 5️⃣ Step 5: 파이프라인 실행

### 설정 확인

```bash
# pfNextflow 디렉토리로 이동
cd pfNextflow

# configs/params.yaml이 있는지 확인
ls -la configs/params.yaml
```

### Docker로 실행 (기본값)

```bash
# 방법 1️⃣: 프로필 지정 없음 (Docker 기본)
nextflow run main.nf

# 방법 2️⃣: 명시적으로 Docker 지정
nextflow run main.nf -profile docker

# 방법 3️⃣: 상세 로그 출력
nextflow run main.nf -profile docker -v

# 방법 4️⃣: 임시 파라미터 오버라이드
nextflow run main.nf -profile docker \
  --do_alignment true \
  --do_ivar false
```

### 실행 중 로그 확인

```bash
# 별도의 터미널에서
tail -f .nextflow.log
```

### 이전 실행 재개

```bash
# 중단된 작업 계속하기
nextflow run main.nf -resume
```

---

## ⚙️ configs/params.yaml 빠른 설정

### BAM 파일 분석 (권장)

```yaml
paths:
  bam_root: "../data/result/03_pf_mapping/mapped_pf"
  bam_glob: "*.bam"
  outdir: "../data/result/analysis_output"

pipeline:
  do_alignment: false    # ✅ BAM이 이미 있음
  do_depth: true
  do_ivar: true
  do_variants: true
  do_afdp: true
```

### FASTQ 분석 (FASTQ → BAM → 분석)

```yaml
paths:
  fastq_root: "../data"
  fastq_glob: "barcode*/*.fastq.gz"
  outdir: "../data/result/analysis_output"

pipeline:
  do_alignment: true     # ✅ FASTQ → BAM 정렬
  do_depth: true
  do_ivar: true
  do_variants: true
```

---

## 🐛 문제 해결

### ❌ "nextflow: command not found"

```bash
# WSL2/Git Bash에서:
curl -s https://get.nextflow.io | bash
chmod +x nextflow
./nextflow -version

# 또는 경로에 추가:
export PATH=$PATH:$(pwd)
nextflow -version
```

### ❌ "Cannot connect to Docker daemon"

```bash
# Docker Desktop 앱이 실행 중인지 확인

# WSL2에서:
sudo service docker start

# 또는 Docker Desktop을 재시작
```

### ❌ "Docker image not found"

```bash
# 이미지가 로드되었는지 확인
docker images | grep biocontainers

# 없으면 다시 로드
docker load -i docker_images/samtools.tar
```

### ❌ "java: command not found"

```powershell
# Java 설치 확인
java -version

# 설치 안 되었으면
choco install openjdk11 -y

# 또는 시스템 재부팅 필요
```

---

## 📝 Windows PowerShell 유용한 명령들

```powershell
# Nextflow 최신 버전 확인
nextflow -version

# 파이프라인 드라이 런 (테스트)
nextflow run main.nf -preview

# 최근 작업 상태 확인
nextflow log

# 작업 캐시 초기화
nextflow clean -f
```

---

## ✅ 설치 완료 확인 체크리스트

- [ ] Java 11+ 설치됨
- [ ] Docker 설치 및 실행 중
- [ ] Nextflow 설치됨
- [ ] Docker 이미지 다운로드/로드됨
- [ ] configs/params.yaml 설정 완료
- [ ] `nextflow run main.nf` 실행 가능

---

## 📚 유용한 리소스

- [Nextflow 공식 문서](https://www.nextflow.io/docs/latest/)
- [Docker 공식 문서](https://docs.docker.com/)
- [WSL2 설정 가이드](https://docs.microsoft.com/en-us/windows/wsl/install)
- 프로젝트 가이드: [REFACTORED_USAGE.md](REFACTORED_USAGE.md)

