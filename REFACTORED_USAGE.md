# 리팩토링된 pfNextflow 사용 가이드

## 📋 변경 사항 요약

이 파이프라인이 **YAML 기반 설정 + Docker 오프라인 실행**으로 리팩토링되었습니다.

### 이전 방식 (❌ 더 이상 권장하지 않음)
```bash
# 매번 main.nf를 수정하고 Nextflow를 설치해야 했음
nextflow run main.nf --bam_root "./data" --bam_glob "*.bam"
```

### 새로운 방식 (✅ 권장)
```bash
# ✅ Nextflow 설치 필요 (→ NEXTFLOW_INSTALL.md 참고)
# ✅ configs/params.yaml 수정만 하면 됨
# ✅ Docker 오프라인 모드로 자동 실행
nextflow run main.nf
```

---

## 🚀 빠른 시작 (5분)

### 1단계: Nextflow 설치
**⚠️ 아직 설치 안 했으면 먼저 하세요!**

[→ NEXTFLOW_INSTALL.md 참고](NEXTFLOW_INSTALL.md)

```bash
# 빠른 체크
java -version      # Java 11+ 필요
docker --version   # Docker 필요
nextflow -version  # Nextflow 필요
```

### 2단계: 기본 구조 확인
```
pfNextflow/
├── configs/
│   └── params.yaml          ← 🔧 이 파일만 수정하세요!
├── data/
│   ├── barcode01/           ← 입력 데이터
│   └── result/
│       ├── 03_pf_mapping/   ← BAM 파일
│       └── analysis_output/ ← 출력 결과
├── main.nf                  ← 메인 파이프라인
└── nextflow.config          ← Docker 기본 설정 (수정 필요 없음)
```

### 3단계: params.yaml 수정

`configs/params.yaml` 파일을 열어서 경로만 확인하세요:

```yaml
# 📁 DATA PATHS
paths:
  bam_root: "../data/result/03_pf_mapping/mapped_pf"
  bam_glob: "*.bam"
  outdir: "../data/result/analysis_output"
```

### 4단계: 시작!

```bash
# ✅ Docker 기본 설정으로 실행
nextflow run main.nf

# ✅ 자동으로 Docker 이미지 다운로드 및 실행
# ✅ 모든 분석이 컨테이너 내에서 수행됨
```

---

## 📝 설정 파일: configs/params.yaml

### 파일 위치
```
pfNextflow/
└── configs/
    └── params.yaml          ← 현재 위치
```

### 주요 설정 섹션

#### 📁 데이터 경로 (paths)
```yaml
paths:
  bam_root: "../data/result/03_pf_mapping/mapped_pf"
  bam_glob: "*.bam"
  fastq_root: "../data"
  fastq_glob: "barcode*/*.fastq.gz"
  outdir: "../data/result/analysis_output"
```

#### 🔬 분석 활성화/비활성화 (pipeline)
```yaml
pipeline:
  do_fastqc: false          # FastQ 품질 검사
  do_alignment: false       # Minimap2 정렬
  do_depth: true            # ✅ 깊이 분석
  do_ivar: true             # ✅ iVar consensus
  do_gene_bins: true        # ✅ 유전자별 분석
  do_variants: true         # ✅ BCFtools 변이 호출
  do_multipanel: false      # 다중 패널 플롯
  do_gene_1x8: false        # 1x8 유전자 플롯
  do_afdp: true             # ✅ AF/DP 피규어
```

#### ⚙️ 파라미터 조정
```yaml
ivar:
  threshold: 0.6            # Consensus 스레시홀드
  min_cons_cov: 10          # 최소 커버리지
  min_var_cov: 1
  genes: "CRT DHFR DHPS K13 MDR1 MSP2 PMI PMIII"

minimap2:
  preset: "sr"              # short reads
  threads: 4

bcftools:
  snp_only: true
  max_depth: 10000
```

---

## 🚀 실행 방법

### Docker 오프라인 모드 (기본값 ✅)

```bash
# 방법 1️⃣: 간단한 실행
nextflow run main.nf

# 방법 2️⃣: 명시적으로 docker 지정
nextflow run main.nf -profile docker

# 방법 3️⃣: 상세 로그 출력
nextflow run main.nf -profile docker -v

# 방법 4️⃣: 중단된 작업 계속
nextflow run main.nf -resume
```

### 커스텀 파라미터로 실행

```bash
# 명령줄에서 파라미터 오버라이드
nextflow run main.nf \
  --do_alignment true \
  --do_ivar false \
  --ivar_threshold 0.8
```

---

## 🎯 일반적인 사용 사례

### Case 1️⃣: 기존 BAM 파일 분석 (가장 일반적)

```yaml
# ✅ 이미 정렬된 BAM 파일이 있는 경우

paths:
  bam_root: "../data/result/03_pf_mapping/mapped_pf"
  bam_glob: "*.bam"
  outdir: "../data/result/analysis_output"

pipeline:
  do_alignment: false       # ✅ BAM이 있으므로 정렬 스킵
  do_depth: true
  do_ivar: true
  do_variants: true
  do_afdp: true
```

**실행:**
```bash
nextflow run main.nf
```

---

### Case 2️⃣: FASTQ에서 처음부터 분석

```yaml
# FASTQ 파일로부터 정렬부터 시작

paths:
  fastq_root: "../data"
  fastq_glob: "barcode*/*.fastq.gz"
  outdir: "../data/result/analysis_output"

pipeline:
  do_fastqc: false         # 선택: 품질 검사
  do_alignment: true       # ✅ FASTQ → BAM 정렬
  do_depth: true
  do_ivar: true
  do_variants: true
```

**실행:**
```bash
nextflow run main.nf
```

---

### Case 3️⃣: 특정 유전자만 고속 분석

```yaml
ivar:
  genes: "K13 MDR1"        # 저항성 유전자만 분석
  threshold: 0.8           # 더 엄격한 스레시홀드

pipeline:
  do_gene_bins: true
  do_multipanel: true
  do_gene_1x8: true
  do_variants: true

# 불필요한 것 비활성화
afdp_figures:
  do_afdp: false
```

---

### Case 4️⃣: 민감도 조정 (저빈도 변이)

```yaml
bcftools:
  snp_only: false          # Indel도 포함
  max_depth: 50000         # 깊이 제한 증가

afdp_figures:
  af_threshold: 0.05       # 낮은 AF에 민감
  min_depth: 5             # 낮은 깊이 허용
```

---

## 📊 출력 구조

모든 분석 결과는 `outdir`에 저장됩니다:

```
analysis_output/
├── depth/                 ← 깊이 분석
├── consensus/             ← iVar consensus 시퀀스
├── gene_bins/             ← 유전자별 깊이
├── calling/
│   └── bcftools/          ← VCF 변이 파일
├── plots/
│   ├── af_dp/             ← AF/DP 피규어
│   ├── gene_bins/         ← 다중 패널 플롯
│   └── gene_1x8/          ← 1x8 유전자 플롯
└── pipeline_info/         ← 파이프라인 리포트
    ├── timeline.html
    └── report.html
```

---

## 📝 팀 & 트릭

### 여러 분석 구성 관리

```bash
# 설정 백업
cp configs/params.yaml configs/params.backup.yaml

# 다른 설정으로 실행
cp configs/params.analysis1.yaml configs/params.yaml
nextflow run main.nf

# 원래대로 복구
cp configs/params.backup.yaml configs/params.yaml
```

### 드라이 런 (실제 실행 없이 테스트)

```bash
nextflow run main.nf -preview
```

### 로그 파일 확인

```bash
# 실행 중 로그 보기
tail -f .nextflow.log

# 이전 실행 로그 보기
nextflow log
```

### 작업 캐시 초기화 (필요시)

```bash
# 모든 캐시 삭제
nextflow clean -f

# 또는 특정 작업만 재실행
nextflow run main.nf -resume
```

---

## 🔧 Docker 오프라인 모드

이 파이프라인은 **모든 작업을 Docker 컨테이너에서 실행**합니다:

### 필요한 Docker 이미지

```
biocontainers/samtools:latest
biocontainers/fastqc:v0.12.1-1-deb_cv1
biocontainers/minimap2:v2.24-1-deb_cv1
andersenlabapps/ivar:latest
biocontainers/bcftools:v1.17-1-deb_cv1
biocontainers/bedtools:v2.31.0-1-deb_cv1
rocker/r-base:4.3.1
```

### 오프라인 환경 설정

[→ OFFLINE_DOCKER.md](OFFLINE_DOCKER.md) 또는 [→ NEXTFLOW_INSTALL.md](NEXTFLOW_INSTALL.md) 참고

```bash
# 이미지 미리 다운로드
docker pull biocontainers/samtools:latest
docker pull biocontainers/fastqc:v0.12.1-1-deb_cv1
# ... (다른 이미지들도)

# 또는 오프라인으로 로드
docker load -i docker_images/samtools.tar
docker load -i docker_images/fastqc.tar
# ... (다른 이미지들도)
```

---

## 📋 파라미터 완전 문서

### paths - 데이터 경로
| 파라미터 | 설명 | 기본값 |
|----------|------|--------|
| `bam_root` | BAM 파일 디렉토리 | "../data/result/03_pf_mapping/mapped_pf" |
| `bam_glob` | BAM 파일 패턴 | "*.bam" |
| `fastq_root` | FASTQ 파일 디렉토리 | "../data" |
| `fastq_glob` | FASTQ 파일 패턴 | "barcode*/*.fastq.gz" |
| `outdir` | 출력 디렉토리 | "../data/result/analysis_output" |

### references - 참조 파일
| 파라미터 | 설명 | 기본값 |
|----------|------|--------|
| `bed` | 증폭자 BED 파일 | amplicons.bed |
| `ref` | 참조 게놈 | PlasmoDB-67_Pfalciparum3D7_Genome.fasta |

### pipeline - 분석 단계
| 파라미터 | 설명 | 기본값 |
|----------|------|--------|
| `do_fastqc` | FastQ 품질 검사 | false |
| `do_alignment` | Minimap2 정렬 | false |
| `do_depth` | 깊이 분석 | true |
| `do_ivar` | iVar consensus | true |
| `do_gene_bins` | 유전자별 분석 | true |
| `do_variants` | BCFtools 변이 | true |
| `do_multipanel` | 다중 패널 플롯 | false |
| `do_gene_1x8` | 1x8 유전자 플롯 | false |
| `do_afdp` | AF/DP 피규어 | true |

### ivar - Consensus 설정
| 파라미터 | 설명 | 기본값 |
|----------|------|--------|
| `threshold` | Consensus 스레시홀드 | 0.6 |
| `min_cons_cov` | 최소 커버리지 | 10 |
| `min_var_cov` | 최소 변이 커버리지 | 1 |
| `genes` | 분석 유전자 | "CRT DHFR DHPS K13 MDR1 MSP2 PMI PMIII" |

### minimap2 - 정렬 설정
| 파라미터 | 설명 | 기본값 |
|----------|------|--------|
| `preset` | 모드 (sr/map-pb/map-ont) | "sr" |
| `threads` | 스레드 수 | 4 |

### bcftools - 변이 호출 설정
| 파라미터 | 설명 | 기본값 |
|----------|------|--------|
| `mode` | 모드 (m/v/mv) | "mv" |
| `snp_only` | SNP만 추출 | true |
| `threads` | 스레드 수 | 4 |
| `max_depth` | 최대 깊이 | 10000 |

### afdp_figures - AF/DP 피규어 설정
| 파라미터 | 설명 | 기본값 |
|----------|------|--------|
| `min_depth` | 최소 깊이 필터 | 15 |
| `af_threshold` | AF 스레시홀드 | 0.50 |
| `shade_xmin` | 음영 X최소 | 0.48 |
| `shade_xmax` | 음영 X최대 | 0.52 |

---

## 🐛 주요 문제 해결

### Nextflow 설치 필요
[→ NEXTFLOW_INSTALL.md](NEXTFLOW_INSTALL.md) 참고

```bash
# 확인
java -version
docker --version
nextflow -version
```

### "No BAM/FASTQ files found" 오류

```yaml
# ❌ 잘못된 경로
paths:
  bam_root: "../data/result/mapping"

# ✅ 올바른 경로 확인
# ls ../data/result/03_pf_mapping/mapped_pf/
# *.bam 파일이 있는지 확인
```

### Docker 이미지 없음

```bash
# 필요한 이미지 다운로드
docker pull biocontainers/samtools:latest
# ... (다른 이미지들)

# 또는 오프라인 로드
docker load -i docker_images/samtools.tar
```

---

## 📞 추가 고급 설정

### 환경별 프로필

```bash
# 로컬 실행 (도구 직접 설치)
nextflow run main.nf -profile local

# SLURM 클러스터 (HPC)
nextflow run main.nf -profile slurm

# Docker (기본값, 권장)
nextflow run main.nf -profile docker
```

### 성능 튜닝

```yaml
resources:
  max_cpus: 32
  max_memory: "128.GB"
  max_time: "240.h"
```

---

## 📚 추가 리소스

- [Nextflow 공식 문서](https://www.nextflow.io/docs/latest/)
- [Nextflow 설치 가이드](NEXTFLOW_INSTALL.md)
- [Docker 오프라인 설정](OFFLINE_DOCKER.md)
- [프로젝트 README](README.md)


---

## 🎯 일반적인 사용 사례

### Case 1: 기존 BAM 파일 분석

```yaml
paths:
  bam_root: "../data/result/03_pf_mapping/mapped_pf"
  bam_glob: "*.bam"
  outdir: "../data/result/analysis_output"

pipeline:
  do_alignment: false    # ✅ BAM 파일이 있으므로 정렬 스킵
  do_depth: true
  do_ivar: true
  do_variants: true
  do_afdp: true
```

**실행:**
```bash
nextflow run main.nf
```

---

### Case 2: FASTQ에서 처음부터 분석

```yaml
paths:
  fastq_root: "../data"
  fastq_glob: "barcode*/*.fastq.gz"
  outdir: "../data/result/analysis_output"

pipeline:
  do_fastqc: false       # 품질 검사 (선택)
  do_alignment: true     # ✅ FASTQ → BAM 정렬
  do_depth: true
  do_ivar: true
  do_variants: true
  do_afdp: true

minimap2:
  preset: "sr"           # 단문제(short reads)
  threads: 4
```

**실행:**
```bash
nextflow run main.nf
```

---

### Case 3: 특정 유전자만 분석

```yaml
ivar:
  genes: "K13 MDR1"      # CRT, DHFR 등의 저항성 유전자만
  
pipeline:
  do_gene_bins: true     # 유전자별 깊이 분석
  do_multipanel: true    # 다중 패널 플롯
  do_gene_1x8: true      # 1x8 유전자 플롯
```

---

### Case 4: 민감도 조정

```yaml
ivar:
  threshold: 0.8         # 더 높은 스레시홀드 (더 엄격)
  min_cons_cov: 20       # 더 높은 커버리지 요구
  min_var_cov: 2         # 더 높은 변이 커버리지

bcftools:
  snp_only: true         # SNP만 추출 (indel 제외)
  max_depth: 5000        # 깊이 제한 변경

afdp_figures:
  af_threshold: 0.10     # 낮은 AF에서도 표시
```

---

## 📊 출력 구조

모든 분석 결과는 `outdir`에 저장됩니다:

```
analysis_output/
├── 01_fastqc/           ← FastQ 품질 검사 (if do_fastqc=true)
├── 02_alignment/        ← 정렬된 BAM 파일 (if do_alignment=true)
├── depth_analysis/      ← 깊이 플롯 (if do_depth=true)
├── consensus/           ← iVar consensus sequences (if do_ivar=true)
├── gene_bins/           ← 유전자별 깊이 분석 (if do_gene_bins=true)
├── calling/
│   └── bcftools/        ← VCF 변이 호출 (if do_variants=true)
├── plots/
│   ├── gene_bins/       ← 다중 패널 플롯 (if do_multipanel=true)
│   ├── gene_1x8/        ← 1x8 유전자 플롯 (if do_gene_1x8=true)
│   └── af_dp/           ← AF/DP 피규어 (if do_afdp=true)
└── pipeline_info/       ← 파이프라인 리포트
```

---

## ⚙️ 파라미터 설명

### paths - 데이터 경로
| 파라미터 | 설명 |
|----------|------|
| `bam_root` | BAM 파일이 있는 디렉토리 |
| `bam_glob` | BAM 파일 이름 패턴 (glob) |
| `fastq_root` | FASTQ 파일이 있는 디렉토리 |
| `fastq_glob` | FASTQ 파일 이름 패턴 |
| `outdir` | 출력 디렉토리 |

### references - 참조 파일
| 파라미터 | 설명 | 기본값 |
|----------|------|--------|
| `bed` | 증폭자 BED 파일 | amplicons.bed |
| `ref` | 참조 게놈 FASTA | PlasmoDB-67_Pfalciparum3D7_Genome.fasta |

### pipeline - 분석 단계 토글
| 파라미터 | 설명 | 기본값 |
|----------|------|--------|
| `do_fastqc` | FastQ 품질 관리 | false |
| `do_alignment` | Minimap2 정렬 | false |
| `do_depth` | 깊이 분석 | true |
| `do_ivar` | iVar consensus 생성 | true |
| `do_gene_bins` | 유전자별 깊이 분석 | true |
| `do_variants` | BCFtools 변이 호출 | true |
| `do_multipanel` | 다중 패널 플롯 | false |
| `do_gene_1x8` | 1x8 유전자 플롯 | false |
| `do_afdp` | AF/DP 피규어 | true |

### ivar - Consensus 파라미터
| 파라미터 | 설명 | 기본값 |
|----------|------|--------|
| `threshold` | 합의 스레시홀드 | 0.6 |
| `min_cons_cov` | 최소 consensus 커버리지 | 10 |
| `min_var_cov` | 최소 변이 커버리지 | 1 |
| `genes` | 분석할 유전자 (공백 분리) | CRT DHFR DHPS K13 MDR1 MSP2 PMI PMIII |

### minimap2 - 정렬 파라미터
| 파라미터 | 설명 | 기본값 |
|----------|------|--------|
| `preset` | 정렬 모드 (sr/map-pb/map-ont/asm5) | sr |
| `threads` | 스레드 수 | 4 |

### bcftools - 변이 호출 파라미터
| 파라미터 | 설명 | 기본값 |
|----------|------|--------|
| `mode` | 모드 (m=multiallelic, v=vcf_outputs) | mv |
| `snp_only` | SNP만 추출 | true |
| `threads` | 스레드 수 | 4 |
| `max_depth` | 최대 깊이 제한 | 10000 |

### afdp_figures - AF/DP 피규어 파라미터
| 파라미터 | 설명 | 기본값 |
|----------|------|--------|
| `min_depth` | 최소 깊이 필터 | 15 |
| `af_threshold` | AF 스레시홀드 | 0.50 |
| `shade_xmin` | 음영 X최소값 | 0.48 |
| `shade_xmax` | 음영 X최대값 | 0.52 |

---

## 🐛 문제 해결

### ❌ "YAML 파일을 찾을 수 없습니다" 오류

resources/params.yaml이 없는 경우입니다.
```bash
# 파일이 있는지 확인
ls -la resources/params.yaml

# 없으면 생성
cp resources/params.yaml.example resources/params.yaml
```

### ❌ "No BAM/FASTQ files found" 오류

params.yaml에서 경로를 확인하세요:
```yaml
paths:
  bam_root: "../data/result/03_pf_mapping/mapped_pf"  # 실제 경로인지 확인
  bam_glob: "*.bam"
```

### ❌ Docker 컨테이너 문제

```bash
# Docker 프로필 사용
nextflow run main.nf -profile docker

# 또는 로컬 실행 (도구 직접 설치)
nextflow run main.nf -profile local
```

---

## 📝 팁 & 트릭

### 여러 분석 구성 관리

여러 설정을 저장해 두셨다가 필요할 때 사용하세요:

```bash
# 설정 파일 백업
cp resources/params.yaml resources/params.analysis1.yaml
cp resources/params.yaml resources/params.analysis2.yaml

# 특정 설정으로 실행
cd resources && cp params.analysis1.yaml params.yaml && cd ..
nextflow run main.nf -profile docker
```

### 드라이 런 (실제 실행 없이 테스트)

```bash
nextflow run main.nf -preview
```

### 상세 로그 보기

```bash
nextflow run main.nf -v 2>&1 | tee run.log
```

### 이전 실행 재개

```bash
nextflow run main.nf -resume
```

---

## 📞 추가 지원

- 이 가이드를 참고하세요: [pfNextflow README.md](README.md)
- Docker 설정: [OFFLINE_DOCKER.md](OFFLINE_DOCKER.md)
