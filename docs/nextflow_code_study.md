# pfNextflow 코드 공부 노트

> `main.nf` 를 중심으로 Nextflow 개념 정리

---

## 1. 설계 철학 — Python과 비교

Python의 `main() + 모듈 패키지` 구조와 사실 동일하다.

| Python | Nextflow |
|--------|----------|
| `def my_func(a, b):` | `process MY_PROCESS { input: ... }` |
| `import utils` | `include { MY_PROCESS } from './modules/my.nf'` |
| `def main():` | `workflow { }` |
| 함수 반환값 | `emit: result_name` |
| 변수에 할당 | `result_ch = MY_PROCESS(input_ch)` |

### 권장 폴더 구조

```
main.nf          ← workflow {} 만 있음  (Python의 main.py)
modules/
  align.nf       ← process ALIGN_READS {}  (Python의 utils/align.py)
  call.nf        ← process CALL_VARIANTS {}
bin/
  script.R       ← 실제 계산 로직  (Python의 utils/helpers.py)
configs/
  params.yaml    ← 설정값 분리
```

### 설계 원칙 3가지

1. **process는 셸 명령 하나만** — `minimap2`, `bcftools call` 각각 분리
2. **파일은 항상 channel로** — 직접 path 문자열 넘기지 말고 `channel.fromPath()`
3. **데이터는 흘려보내기** — push하지 말고 channel 연산자(`map`, `join`, `groupTuple`)로 변환

---

## 2. workflow vs process — 역할 분리

```nextflow
// process = 일꾼. 입력받아 셸 실행, 결과 emit
process CALL_VARIANTS {
    input:  tuple val(barcode), path(bam)
    output: path "*.vcf.gz", emit: vcf
    script:
    """
    bcftools mpileup ... | bcftools call ...
    """
}

// workflow = 지휘자. process 조립, channel 연결
workflow {
    bam_ch = channel.fromPath("data/**/*.bam")
    vcf_ch = CALL_VARIANTS(bam_ch)   // 반환값 필요하면 변수 할당
    PF_DEPTH(bam_ch, ...)             // 결과가 publishDir만 가면 변수 불필요
}
```

**변수에 할당하는 경우 vs 안 하는 경우:**

```nextflow
vcf_ch = CALL_VARIANTS(...)   // 아래에서 vcf_ch 를 또 쓸 때
PF_DEPTH(bam_ch, ...)         // 결과가 publishDir 로만 가고 끝날 때
```

---

## 3. Channel 핵심 연산자

### `channel.value()` vs `channel.fromPath()`

```nextflow
// value channel: 무한 재사용 (공통 리소스)
ref_ch = channel.value(file("resources/ref.fasta"))
// → 모든 샘플에 같은 ref 파일을 재사용 가능

// queue channel: 1번만 소비됨 (샘플별 파일)
bam_ch = channel.fromPath("data/**/*.bam")
// → 각 BAM이 한 번씩 흘러감
```

### `.map { }` — 형태 변환

```nextflow
// Path 객체를 tuple 로 변환
channel.fromPath("data/**/*.bam")
  .map { bam -> tuple(bam.name.replaceAll(/\..+$/, ""), bam) }

// 결과:
// ["barcode01", Path(barcode01.sorted.bam)]
// ["barcode02", Path(barcode02.sorted.bam)]
```

> `{ bam -> ... }` 에서 `bam` 은 임의 변수명 (Python lambda의 파라미터와 동일)

### `.groupTuple()` — 같은 key끼리 묶기

```nextflow
// 입력:
// ["barcode01", file_a]
// ["barcode01", file_b]
// ["barcode02", file_c]

.groupTuple()

// 출력:
// ["barcode01", [file_a, file_b]]
// ["barcode02", [file_c]]
```

### `.flatMap()` vs `.flatten()`

```nextflow
// flatMap: 아이템마다 여러 개로 펼치기
vcf_ch
  .flatMap { dir -> file("${dir}/**/*.vcf.gz") }
// 디렉토리 1개 → 그 안의 vcf 파일 여러 개로 확장

// flatten: 리스트 아이템을 개별 아이템으로 분리
SPLICE_TRANSLATE.out.nt_fastas  // 샘플별로 [파일, 파일, ...] 리스트
  .flatten()                    // 파일 1개씩으로 분리
```

Python 비유:
```python
# flatMap ≈
[f for dir in dirs for f in glob(dir + "/**/*.vcf.gz")]

# flatten ≈
itertools.chain.from_iterable(list_of_lists)
```

### `.join()` — 두 채널을 key로 합치기

```nextflow
nt_by_gene_ch  // ["CRT", [nt파일들...]]
  .join(aa_by_gene_ch)  // ["CRT", [aa파일들...]]

// 결과:
// ["CRT", [nt파일들...], [aa파일들...]]
```

SQL의 `INNER JOIN ON key` 와 동일.

### `.collect()` — 전체 모아서 한 번에

```nextflow
all_vcfs_ch = vcf_ch.collect()
// 모든 vcf가 다 나올 때까지 기다렸다가 한 번에 다음 process로
```

### `.ifEmpty { error }` — 빈 채널 방어

```nextflow
channel.fromPath(...)
  .ifEmpty { error "❌ 파일이 없습니다" }
```

---

## 4. 실시간 출력 — log.info vs process echo

### 문제

```nextflow
// 이건 workflow 정의 시점에 출력됨 (실행 전에 다 찍힘)
log.info "Running: CDS Splice, Translate & Variant Table"
```

### 해결책

**① process script 안에 echo (가장 실용적)**

```nextflow
process SPLICE_TRANSLATE {
    script:
    """
    echo "▶ [SPLICE_TRANSLATE] 시작: ${barcode}"
    python3 splice_and_translate.py ...
    echo "✅ [SPLICE_TRANSLATE] 완료: ${barcode}"
    """
}
```

**② nextflow.config 에 `echo = true`**

```groovy
process {
    echo = true   // process의 stdout을 터미널에 출력
}
```

**③ `.view()` — channel 완료 감지**

```nextflow
vcf_ch = CALL_VARIANTS(...)
vcf_ch.view { dir -> "✅ CALL_VARIANTS 완료: ${dir.name}" }
// process output이 emit될 때 = task 완료 시점에 출력
```

**④ trace/timeline report (사후 분석)**

```groovy
// nextflow.config
timeline { enabled = true }
trace    { enabled = true }
report   { enabled = true }
```

| 방법 | 타이밍 | 용도 |
|------|--------|------|
| `log.info` | 정의 시점 (항상 앞에) | 설정값 확인 |
| script `echo` | 실제 실행 시점 | 실시간 진행 |
| `.view()` | task 완료 시점 | 완료 감지 |
| trace/report | 사후 | 성능 분석 |

---

## 5. SPLICE → ALIGN 블록 전체 흐름 해설

```nextflow
// ── 입력: 샘플별 디렉토리 ──────────────────────────────────────────
vcf_ch
  .map { dir -> tuple(dir.name, dir) }
  .set { variants_with_barcode_ch }
// ["barcode01", Path(barcode01/)]
// ["barcode02", Path(barcode02/)]

// ── SPLICE_TRANSLATE: 샘플별 실행 ─────────────────────────────────
SPLICE_TRANSLATE(variants_with_barcode_ch, ...)
// 출력 (샘플별로 리스트):
// out.nt_fastas: [barcode01.CRT.nt.fasta, barcode01.DHFR.nt.fasta, ...]
// out.nt_fastas: [barcode02.CRT.nt.fasta, barcode02.DHFR.nt.fasta, ...]

// ── flatten: 리스트 → 개별 파일 ────────────────────────────────────
SPLICE_TRANSLATE.out.nt_fastas
  .flatten()
// barcode01.CRT.nt.fasta
// barcode01.DHFR.nt.fasta
// barcode02.CRT.nt.fasta
// ...

// ── 유전자별 재그룹핑 ──────────────────────────────────────────────
  .filter { f -> !f.name.startsWith('EMPTY') }  // 빈 파일 제거
  .unique  { it.name }                           // -resume 시 중복 방지
  .map     { f -> tuple(f.name.tokenize('.')[1], f) }
  // "barcode01.CRT.nt.fasta".tokenize('.') = ["barcode01", "CRT", "nt", "fasta"]
  // [1] = "CRT"
  // → ["CRT", barcode01.CRT.nt.fasta]
  .groupTuple()
  // → ["CRT",  [barcode01.CRT.nt.fasta, barcode02.CRT.nt.fasta, ...]]
  // → ["DHFR", [barcode01.DHFR.nt.fasta, barcode02.DHFR.nt.fasta, ...]]
  .set { nt_by_gene_ch }

// ── NT + AA join ────────────────────────────────────────────────────
nt_by_gene_ch
  .join(aa_by_gene_ch)
  // → ["CRT",  [nt파일들...], [aa파일들...]]
  .set { by_gene_ch }

// ── ALIGN_VARIANTS: 유전자별 실행 ──────────────────────────────────
ALIGN_VARIANTS(by_gene_ch, ...)
```

### 변환 흐름 그림

```
샘플 단위                   flatten 후                 유전자 단위
──────────────────          ──────────────────          ──────────────────────────
barcode01/ ─SPLICE─► [CRT.nt, DHFR.nt]                CRT.nt  ─map─► ("CRT", f)
barcode02/ ─SPLICE─► [CRT.nt, DHFR.nt]  ─flatten─►   DHFR.nt ─map─► ("DHFR", f)
                                                        ...
                                                              ─groupTuple─► ("CRT", [모든샘플])
                                                                            ("DHFR", [모든샘플])
                                                                    ─ALIGN_VARIANTS─►
```

**핵심:** SPLICE는 샘플 단위, ALIGN은 유전자 단위로 돌기 때문에 그 사이에 `flatten → groupTuple` 로 재편성이 필요하다.

---

## 6. 디버깅 흐름

```
에러 발생
  ↓
.nextflow.log 에서 실패 task hash 확인
  ↓
cd work/<2자리>/<나머지해시>/
  ↓
cat .command.err        # 에러 메시지
bash .command.sh        # 셸 스크립트 직접 재실행 (환경 그대로)
  ↓
문제 수정
  ↓
nextflow run main.nf -resume   # 성공한 task는 캐시 재사용
```

### work/ 폴더 구조 이해

```
실행 1회
  └── process ALIGN_READS
        ├── barcode01 → work/3a/f7b2c0.../   ← task 1개
        ├── barcode02 → work/1c/a8d391.../   ← task 1개
  └── process CALL_VARIANTS
        ├── barcode01 → work/7f/cc1234.../
        ...
```

- 폴더 이름 = `입력파일 + 스크립트 내용`의 SHA256 hash
- 같은 입력이면 같은 hash → `-resume` 캐시 재사용 원리
- **로그 1개 ≠ work 폴더 1개** : 24샘플 × 10 process = ~240개 폴더

### 어느 폴더가 어느 task?

```bash
nextflow log                              # 최근 실행 목록
nextflow log <run_name> -f name,workdir  # task별 work 경로 출력
grep "barcode01" .nextflow.log | grep "work/"  # 로그에서 직접 검색
```

### `.nextflow.log` 가 계속 쌓이는 이유

실행마다 새 `.nextflow.log` 를 만들고 기존 건 `.nextflow.log.1`, `.2` 로 rotate (최대 10개 유지). 이전 실행 기록 비교용으로 의도적인 설계.

---

## 7. 유용한 디버깅 도구 모음

```nextflow
// ① channel 내용 확인
my_ch.view { "DEBUG: $it" }

// ② dry run (실제 실행 없이 연결만 검증)
// nextflow run main.nf -stub

// ③ 상세 로그
// NXF_DEBUG=1 nextflow run main.nf
```

```bash
# ④ 특정 process만 재실행
nextflow run main.nf -resume -entry MY_PROCESS

# ⑤ trace 파일로 실행 시간 분석
# results/pipeline_info/trace.txt 참고
```
