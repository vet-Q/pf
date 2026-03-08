# 🚀 pfNextflow - 빠른 시작 (Windows)

## 📋 파일 설명

### 1️⃣ 필수 실행 파일

| 파일 | 목적 | 언제 | 비고 |
|------|------|------|------|
| **install_nextflow.bat** | Nextflow 설치 | 첫 1회만 | 관리자 권한 필요 |
| **run_pipeline.bat** | 파이프라인 실행 | 매번 | 더블클릭으로 실행 |

### 2️⃣ 선택 실행 파일

| 파일 | 목적 | 언제 |
|------|------|------|
| **download_docker_images.bat** | Docker 이미지 미리 다운로드 | 첫 실행 전 (선택) |
| **setup.bat** | 종합 설정 확인 | 문제 발생 시 |

---

## 🎯 빠른 시작 가이드 (5분)

### 📍 Step 1: Nextflow 설치 (첫 1회만)

1. **현재 폴더를 탐색기에서 열기**
   ```
   C:\Users\kwono\OneDrive\문서\python project\pfNextflow
   ```

2. **`install_nextflow.bat` 더블클릭**
   - Java 확인
   - Nextflow 다운로드 및 설치
   - 완료 시 자동으로 닫힘

3. **확인 메시지가 나올 때까지 기다리기**

### 📍 Step 2: Docker 이미지 준비 (선택사항)

필요한 Docker 이미지를 미리 다운로드하려면:

1. **`download_docker_images.bat` 더블클릭** (약 5-10분)
   - 자동으로 모든 이미지 다운로드

또는 첫 실행 시에 자동으로 다운로드됩니다.

### 📍 Step 3: 설정 파일 수정 (중요!)

1. **`configs/params.yaml` 파일 열기**
   - 텍스트 에디터로 열기

2. **데이터 경로 확인/수정**
   ```yaml
   paths:
     bam_root: "../data/result/03_pf_mapping/mapped_pf"  # 확인!
     bam_glob: "*.bam"
     outdir: "../data/result/analysis_output"
   ```

3. **필요시 분석 옵션 변경**
   ```yaml
   pipeline:
     do_depth: true       # ✅ 활성화
     do_ivar: true        # ✅ 활성화
     do_variants: true    # ✅ 활성화
     do_afdp: true        # ✅ 활성화
   ```

4. **저장 (Ctrl+S)**

### 📍 Step 4: 파이프라인 실행!

1. **`run_pipeline.bat` 더블클릭**
   - Java 설정 자동 확인
   - Nextflow 파이프라인 시작
   - Docker 컨테이너에서 분석 진행

2. **진행 상황 모니터링**
   - 터미널에 로그가 출력됨
   - 완료 시 결과가 `configs/params.yaml`의 `outdir`에 저장됨

---

## ❓ 대기 화면 (이러한 메시지가 나오면 정상)

설치 또는 실행 중에 "계속하려면 아무 키나 누르세요..." 메시지가 나오면:
- **아무 키나 누르기** (예: Enter)

---

## 🐛 문제 해결

### ❌ "`java` 명령을 찾을 수 없습니다"

**해결:**
1. 모든 CMD/PowerShell 창 닫기
2. 윈도우 **시스템 재부팅**
3. 다시 `install_nextflow.bat` 실행

### ❌ "Docker를 찾을 수 없습니다"

**해결:**
1. Docker Desktop 설치 확인: https://www.docker.com/products/docker-desktop
2. Docker Desktop 앱 실행 (백그라운드에서 대기)
3. Docker가 실행 중인지 확인: 
   - 시스템 트레이에서 Docker 아이콘 확인

### ❌ "configs/params.yaml를 찾을 수 없습니다"

**해결:**
1. 현재 폴더가 `pfNextflow`인지 확인
2. 폴더 안에 `configs` 폴더가 있는지 확인
3. `configs` 폴더 안에 `params.yaml`이 있는지 확인

### ❌ "파이프라인 오류 발생"

**로그 확인:**
1. 같은 폴더에 `.nextflow.log` 파일 확인
2. 에러 메시지 읽고 경로/설정 수정
3. `run_pipeline.bat` 다시 실행

---

## 📚 고급 설명

### 🔍 설정 파일 상세 정보

`configs/params.yaml`에서 조정 가능한 주요 파라미터:

```yaml
# 입력/출력 경로
paths:
  bam_root: "../data/result/03_pf_mapping/mapped_pf"  # BAM 파일 디렉토리
  bam_glob: "*.bam"                                    # 파일 이름 패턴
  outdir: "../data/result/analysis_output"             # 결과 저장 위치

# iVar 파라미터
ivar:
  threshold: 0.6            # Consensus 스레시홀드 (0.0 ~ 1.0)
  min_cons_cov: 10          # 최소 커버리지
  genes: "CRT DHFR DHPS"    # 분석할 유전자 (공백 분리)

# 실행할 분석 선택
pipeline:
  do_depth: true            # 깊이 분석
  do_ivar: true             # Consensus 생성
  do_variants: true         # 변이 호출
  do_afdp: true             # AF/DP 플롯
```

### 📊 결과 파일 위치

`configs/params.yaml`에서 지정한 `outdir` 위치:

```
analysis_output/
├── depth/                  ← 깊이 분석
├── consensus/              ← Consensus 시퀀스
├── calling/
│   └── bcftools/           ← VCF 파일
├── plots/
│   ├── af_dp/              ← AF/DP 플롯
│   └── gene_bins/          ← 유전자별 플롯
└── pipeline_info/          ← 파이프라인 리포트
```

### 🔄 재실행 및 이어서 실행

같은 폴더에서 다시 실행:
```
run_pipeline.bat  # 처음부터 실행
```

중단된 작업 계속하기:
1. 텍스트 에디터에서 `nextflow run main.nf -resume` 실행
2. 또는 PowerShell/CMD에 직접 입력

---

## 📞 추가 정보

- [Nextflow 공식 문서](https://www.nextflow.io/docs/latest/)
- [Docker 설치 가이드](https://docs.docker.com/desktop/install/windows-install/)
- 프로젝트 상세 가이드: REFACTORED_USAGE.md

---

## ✅ 체크리스트

- [ ] Java 설치됨 (`java -version` 실행 가능)
- [ ] Docker Desktop 설치하고 실행 중
- [ ] `install_nextflow.bat` 실행 완료
- [ ] `configs/params.yaml` 경로 확인/수정됨
- [ ] `run_pipeline.bat` 준비 완료

**모두 완료되었으면 `run_pipeline.bat` 더블클릭하세요! 🚀**
