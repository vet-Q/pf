# pfNextflow 실행 가이드

## 준비

### Docker 이미지 준비 (오프라인 환경)

`containers/` 폴더에 다음 이미지를 다운받아 넣으세요:
```
containers/fastqc.tar.gz
containers/minimap2.tar.gz
containers/bcftools.tar.gz
containers/samtools.tar.gz
containers/bedtools.tar.gz
containers/ivar.tar.gz
containers/r-base.tar.gz
```

이미지 로드:
```bash
cd containers/
for file in *.tar.gz; do
    docker load < "$file"
done
```

## 실행

### 기본 실행 (FastQC + 정렬 + 변이)
```bash
nextflow run main.nf \
  --fastq_root D:\input \
  --outdir D:\output \
  --do_fastqc true \
  --do_alignment true \
  --do_variants true \
  -profile docker
```

### 추가 옵션
```bash
# AF/DP 그래프
--do_afdp true

# 모든 분석
--do_fastqc true \
--do_alignment true \
--do_variants true \
--do_afdp true \
--do_gene_bins true \
--do_multipanel true \
--do_gene_1x8 true

# 계속하기 (실패한 지점부터)
-resume
```

## 결과
```
D:\output\
├── qc/fastqc/        품질 보고서
├── alignment/        정렬 파일
├── calling/          VCF 변이
└── pipeline_info/    실행 로그
```

## 정리
```bash
# 캐시 제거
nextflow clean -f
```
