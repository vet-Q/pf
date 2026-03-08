# Docker 이미지 직접 사용 (오프라인)

## 폴더 구조
```
pfNextflow/
└── containers/
    ├── fastqc.tar.gz
    ├── minimap2.tar.gz
    ├── bcftools.tar.gz
    ├── samtools.tar.gz
    ├── bedtools.tar.gz
    ├── ivar.tar.gz
    └── r-base.tar.gz
```

## 이미지 로드

```bash
# containers 폴더의 모든 이미지 로드
cd containers/
for file in *.tar.gz; do
    docker load < "$file"
done

# 확인
docker images
```

## Nextflow 실행

```bash
cd pfNextflow
nextflow run main.nf \
  --fastq_root D:\input \
  --outdir D:\output \
  --do_fastqc true \
  --do_alignment true \
  --do_variants true \
  -profile docker
```

완료! 로컬 Docker 이미지를 사용합니다.
