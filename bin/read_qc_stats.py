#!/usr/bin/env python3
"""
read_qc_stats.py  [pysam-free — samtools view 파이프 방식]
===========================================================
Per-amplicon, per-barcode read quality statistics for P. falciparum
Nanopore long-read amplicon sequencing.

pysam 없이 'samtools view' stdout을 직접 파싱합니다.
sankey_efficiency.py와 동일한 방식이므로 WSL/Linux에서 samtools만 있으면 됩니다.

── 계산 지표 4가지 ──────────────────────────────────────────────────────────────
  1. Read length 분포      — 유전자별 stripplot + boxplot (amplicon/CDS 기준선)
  2. ORF full-span rate    — read가 CDS 전체 스팬을 커버하는 비율 (bar chart)
  3. AT 함량 (%)           — Pf AT-rich 특성 확인 (boxplot, 81% 기준선)
  4. Homopolymer 밀도 (%)  — ≥3/6/9 bp 연속 동종염기 비율 (grouped bar)

── 입력 ────────────────────────────────────────────────────────────────────────
  --bam        정렬된 BAM 파일
  --bed        amplicons.bed  (chrom  start  end  gene_name, 탭 구분, 0-based)
  --cds        cds_coords.tsv (gene  chrom  cds_start  cds_end  strand  exon_rank)
  --out-prefix 출력 파일 공통 접두사
  --min-mapq   최소 MAPQ (기본: 0)
  --min-qlen   최소 read 길이 bp (기본: 500)

── 출력 ────────────────────────────────────────────────────────────────────────
  {prefix}_stats.tsv   — read별 통계 테이블
  {prefix}.pdf         — 4-panel matplotlib 그림

── 의존성 ──────────────────────────────────────────────────────────────────────
  Python: pandas, matplotlib (pysam 불필요)
  시스템: samtools (PATH에 있어야 함)
  conda install pandas matplotlib numpy

── 사용 예 ────────────────────────────────────────────────────────────────────
  python3 bin/read_qc_stats.py \\
      --bam  results/alignment/minimap2/barcode01.Pf3D7.sorted.bam \\
      --bed  resources/amplicons.bed \\
      --cds  resources/cds_coords.tsv \\
      --out-prefix results/qc/barcode01_read_qc
"""

import argparse
import os
import re
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")   # 화면이 없는 서버(WSL/HPC) 환경에서도 그림 파일 생성 가능
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def _check_samtools():
    """samtools가 PATH에 있는지 확인하고 없으면 명확한 오류로 종료."""
    if shutil.which("samtools") is None:
        sys.exit(
            "[ERROR] samtools를 찾을 수 없습니다.\n"
            "        WSL/Linux 환경에서 실행하거나 samtools를 PATH에 추가하세요.\n"
            "        예시: wsl bash -c \"python3 bin/read_qc_stats.py ...\""
        )


_CIGAR_REF_OPS = re.compile(r'(\d+)([MDNX=])')


def _ref_end_0based(pos_1based, cigar):
    """
    SAM POS(1-based)와 CIGAR 문자열로 0-based 참조 끝 좌표(exclusive)를 계산합니다.
    pysam의 reference_end와 동일한 값을 반환합니다.
    """
    ref_consuming = sum(int(n) for n, op in _CIGAR_REF_OPS.findall(cigar))
    return (pos_1based - 1) + ref_consuming  # 0-based exclusive


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 1. 좌표 파일 로드 함수들
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def load_amplicons(bed_path):
    """
    amplicons.bed를 읽어 유전자별 앰플리콘 좌표 딕셔너리를 반환합니다.

    BED 포맷은 0-based half-open [start, end) 입니다.
    예: Pf3D7_07_v3  402384  406341  CRT

    반환: { "CRT": ("Pf3D7_07_v3", 402384, 406341), ... }
    """
    amplicons = {}
    with open(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            chrom = parts[0]
            start = int(parts[1])   # 0-based 시작 (포함)
            end   = int(parts[2])   # 0-based 끝 (미포함)
            gene  = parts[3]
            amplicons[gene] = (chrom, start, end)
    return amplicons


def load_cds_spans(cds_path):
    """
    cds_coords.tsv에서 유전자별 CDS 전체 게놈 범위를 계산합니다.

    CRT처럼 다중 엑손 유전자는 엑손이 여러 행이므로,
    min(cds_start) ~ max(cds_end)를 전체 CDS 스팬으로 사용합니다.
    (엑손 사이의 인트론 구간은 포함되지만, ORF span 판단에는 이 범위로 충분합니다.)

    0-based 좌표 반환 (BED-style).
    반환: { "CRT": ("Pf3D7_07_v3", 403221, 405951, "+"), ... }
    """
    df = pd.read_csv(cds_path, sep="\t", comment="#")

    # 유전자별로 집계: 첫 번째 chrom, 전체 cds 중 min start, max end, strand
    spans = (
        df.groupby("gene")
          .agg(
              chrom   = ("chrom",     "first"),
              cds_min = ("cds_start", "min"),
              cds_max = ("cds_end",   "max"),
              strand  = ("strand",    "first"),
          )
          .reset_index()
    )

    result = {}
    for _, row in spans.iterrows():
        result[row["gene"]] = (
            row["chrom"],
            int(row["cds_min"]),
            int(row["cds_max"]),
            row["strand"],
        )
    return result


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 2. 시퀀스 품질 계산 함수들
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def compute_at_percent(seq):
    """
    DNA 서열의 AT 함량(%)을 반환합니다.

    Plasmodium falciparum 게놈은 전체적으로 ~81% AT rich이며,
    CDS 영역도 60-80% A+T입니다.
    N 염기는 분모에 포함(총 길이 기준)하여 보수적으로 계산합니다.

    매개변수:
        seq (str): DNA 서열 문자열

    반환:
        float: AT 함량 (%) — 서열이 없으면 NaN
    """
    if not seq:
        return float("nan")
    seq_upper = seq.upper()
    at_count  = seq_upper.count("A") + seq_upper.count("T")
    return 100.0 * at_count / len(seq_upper)


def compute_homopolymer_density(seq, min_len):
    """
    시퀀스 내에서 길이 min_len 이상인 homopolymer run이 차지하는 비율(%)을
    반환합니다.

    Homopolymer run이란 같은 염기가 min_len개 이상 연속되는 구간을 말합니다.
    예: "AAATTTGGGGCC" (12 bp)
        min_len=3: AAA(3) + TTT(3) + GGGG(4) = 10 bp → 83.3%
        min_len=6: 없음 → 0%

    Nanopore 시퀀서의 ONT R9.4 / R10 케미스트리는 동종염기 반복 구간에서
    삽입/삭제 오류가 많이 발생합니다. 특히 POLY-A/T 구간이 긴 Pf 게놈에서
    DHFR intron 및 K13 5'UTR 근처가 문제가 됩니다.

    매개변수:
        seq     (str): DNA 서열
        min_len (int): 최소 homopolymer 길이 기준 (예: 3, 6, 9)

    반환:
        float: homopolymer 비율 (%)
    """
    if not seq:
        return 0.0
    seq_upper        = seq.upper()
    total_hp_bases   = 0
    i = 0
    while i < len(seq_upper):
        j = i + 1
        # 같은 염기가 이어지는 한 j를 전진
        while j < len(seq_upper) and seq_upper[j] == seq_upper[i]:
            j += 1
        run_len = j - i
        if run_len >= min_len:
            total_hp_bases += run_len
        i = j
    return 100.0 * total_hp_bases / len(seq_upper)


def _spans_full_cds(ref_start_0based, ref_end_0based_excl, cds_start, cds_end):
    """
    read의 0-based 참조 좌표가 CDS 전체 범위를 완전히 커버하는지 확인합니다.

    판단 기준:
        ref_start <= cds_start  AND  ref_end_excl >= cds_end

    매개변수:
        ref_start_0based      : read 정렬 시작 (0-based, 포함) = POS - 1
        ref_end_0based_excl   : read 정렬 끝 (0-based, 미포함) = POS-1 + ref_bases
        cds_start             : CDS 시작 (0-based)
        cds_end               : CDS 끝 (0-based, 미포함)

    반환:
        bool
    """
    return ref_start_0based <= cds_start and ref_end_0based_excl >= cds_end


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 3. BAM 파싱 메인 함수
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def collect_read_stats(bam_path, amplicons, cds_spans, min_mapq=0, min_qlen=500):
    """
    BAM 파일을 파싱하여 앰플리콘별 리드 통계를 수집합니다.

    pysam 대신 'samtools view' subprocess를 사용하여 SAM 텍스트를 직접 파싱합니다.
    samtools_efficiency.py와 동일한 방식이므로 WSL/Linux의 samtools만 있으면 됩니다.

    samtools 필터 플래그:
      -F 2304 = secondary(0x100) + supplementary(0x800) 제외
      -q min_mapq = 최소 MAPQ
      region = chrom:start-end (1-based, SAM 표준)

    SAM 컬럼 인덱스 (0-based):
      col[1]  FLAG   (int)
      col[3]  POS    (1-based ref 시작)
      col[5]  CIGAR
      col[9]  SEQ

    매개변수:
        bam_path  : BAM 파일 경로 (str)
        amplicons : { gene: (chrom, start_0based, end_0based) }
        cds_spans : { gene: (chrom, cds_min_0based, cds_max_0based, strand) }
        min_mapq  : 최소 mapping quality (기본 0)
        min_qlen  : 최소 read 길이 bp (기본 500)

    반환:
        pd.DataFrame: gene, read_len, spans_orf, at_pct, hp3_pct, hp6_pct, hp9_pct
    """
    _check_samtools()
    records = []

    for gene, (chrom, amp_start, amp_end) in amplicons.items():

        # ── CDS 전체 범위 결정 ────────────────────────────────────────────────
        if gene in cds_spans:
            _, cds_s, cds_e, _ = cds_spans[gene]
        else:
            print(f"[WARN] {gene}: cds_coords.tsv에 없음. amplicon 범위를 CDS 대신 사용합니다.",
                  file=sys.stderr)
            cds_s, cds_e = amp_start, amp_end

        # ── samtools view 호출 ───────────────────────────────────────────────
        # BED는 0-based이므로 samtools region(1-based)으로 변환: start+1 ~ end
        region = f"{chrom}:{amp_start + 1}-{amp_end}"
        cmd = [
            "samtools", "view",
            "-F", "2304",          # secondary + supplementary 제외
            "-q", str(min_mapq),   # 최소 MAPQ
            bam_path,
            region,
        ]
        try:
            result = subprocess.run(
                cmd, capture_output=True, text=True, check=True
            )
        except subprocess.CalledProcessError as e:
            print(f"[WARN] {gene}: samtools view 실패 — {e.stderr.strip()}",
                  file=sys.stderr)
            continue

        # ── SAM 텍스트 파싱 ──────────────────────────────────────────────────
        n_total = n_pass = 0
        for line in result.stdout.splitlines():
            if not line or line.startswith("@"):
                continue
            n_total += 1
            cols = line.split("\t")
            if len(cols) < 10:
                continue

            seq   = cols[9]
            cigar = cols[5]

            # CIGAR가 '*'이면 맵핑 없음 → skip
            if cigar == "*" or seq == "*":
                continue

            read_len = len(seq)
            if read_len < min_qlen:
                continue

            # POS(1-based) → ref start/end (0-based)
            pos_1based   = int(cols[3])
            ref_start    = pos_1based - 1
            ref_end_excl = _ref_end_0based(pos_1based, cigar)

            # ── 4가지 지표 계산 ───────────────────────────────────────────────
            full_span = _spans_full_cds(ref_start, ref_end_excl, cds_s, cds_e)
            at_pct    = compute_at_percent(seq)
            hp3_pct   = compute_homopolymer_density(seq, min_len=3)
            hp6_pct   = compute_homopolymer_density(seq, min_len=6)
            hp9_pct   = compute_homopolymer_density(seq, min_len=9)

            records.append({
                "gene"      : gene,
                "read_len"  : read_len,
                "spans_orf" : full_span,
                "at_pct"    : at_pct,
                "hp3_pct"   : hp3_pct,
                "hp6_pct"   : hp6_pct,
                "hp9_pct"   : hp9_pct,
            })
            n_pass += 1

        print(f"[INFO]   {gene}: {n_pass}/{n_total} reads passed filters",
              file=sys.stderr)

    return pd.DataFrame(records)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 4. 시각화
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def plot_read_qc(df, amplicons, cds_spans, out_pdf, barcode_label=""):
    """
    4-panel QC 그림을 생성하여 PDF(또는 다른 확장자)로 저장합니다.

    Panel 배치:
        [A: Read length stripplot] [B: ORF full-span rate %]
        [C: AT content % boxplot ] [D: Homopolymer density  ]

    매개변수:
        df            : collect_read_stats() 반환 DataFrame
        amplicons     : load_amplicons() 반환 딕셔너리
        cds_spans     : load_cds_spans() 반환 딕셔너리
        out_pdf       : 출력 파일 경로 (문자열, 확장자에 따라 PNG/PDF 자동 결정)
        barcode_label : 그림 제목에 표시할 barcode 이름 (예: "barcode01")
    """
    # 유전자 목록 — amplicons.bed에 정의된 순서 유지
    genes    = list(amplicons.keys())
    n_genes  = len(genes)
    gene_idx = {g: i for i, g in enumerate(genes)}

    # ── 색상 팔레트 ──────────────────────────────────────────────────────────
    # tab10은 최대 10색. 유전자가 10개 이상이면 순환합니다.
    cmap   = plt.get_cmap("tab10")
    colors = [cmap(i % 10) for i in range(n_genes)]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(
        f"Read QC Statistics — {barcode_label}",
        fontsize=14, fontweight="bold", y=1.01
    )

    ax_len, ax_span = axes[0]
    ax_at,  ax_hp   = axes[1]

    # ══════════════════════════════════════════════════════════════════════════
    # Panel A: Read Length Stripplot
    # ──────────────────────────────────────────────────────────────────────────
    # 각 유전자/앰플리콘에 대해 read 길이를 점으로 표시합니다.
    # 점이 겹치지 않도록 x축 방향으로 jitter(무작위 미세 이동)를 추가합니다.
    # 기준선 두 가지:
    #   - 회색 점선 : amplicon 전체 범위 길이 (BED end - start)
    #   - 삼각형 마커: CDS 스팬 길이 (intron 포함 genomic span)
    # ══════════════════════════════════════════════════════════════════════════
    # seed 고정으로 jitter 재현성 확보
    rng = np.random.default_rng(seed=42)

    for gene in genes:
        sub = df[df["gene"] == gene]
        if sub.empty:
            continue
        xi  = gene_idx[gene]
        # ±0.20 범위의 균일 분포 jitter — 점들이 좌우로 분산되어 밀도 파악 용이
        jit = rng.uniform(-0.20, 0.20, size=len(sub))
        ax_len.scatter(
            xi + jit,
            sub["read_len"],
            s=6,                        # 점 크기 (작을수록 밀집 시 가독성 향상)
            alpha=0.35,                 # 반투명으로 겹침 표현
            color=colors[xi],
            linewidths=0,
            zorder=2,
        )

    # Amplicon 기대 길이 기준선 (회색 점선)
    for i, gene in enumerate(genes):
        _, amp_s, amp_e = amplicons[gene]
        amp_len = amp_e - amp_s
        ax_len.plot(
            [i - 0.35, i + 0.35], [amp_len, amp_len],
            color="dimgray", lw=1.5, ls="--", zorder=3
        )
        # amplicon 길이 텍스트 주석 (점선 바로 위)
        ax_len.text(i, amp_len * 1.01, f"{amp_len:,}", ha="center",
                    fontsize=6, color="dimgray")

    # CDS 스팬 길이 기준선 (색상 삼각형 마커)
    for i, gene in enumerate(genes):
        if gene not in cds_spans:
            continue
        _, cs, ce, _ = cds_spans[gene]
        cds_len = ce - cs
        ax_len.scatter(
            i, cds_len,
            marker="^", s=100, color=colors[i],
            edgecolors="black", linewidths=0.7, zorder=4
        )

    ax_len.set_xticks(range(n_genes))
    ax_len.set_xticklabels(genes, rotation=40, ha="right")
    ax_len.set_ylabel("Read length (bp)")
    ax_len.set_title("A. Read Length Distribution", loc="left", fontweight="bold")
    # 범례
    legend_handles = [
        Line2D([0], [0], color="dimgray", lw=1.5, ls="--", label="Amplicon span"),
        Line2D([0], [0], marker="^", color="gray", lw=0,
               markeredgecolor="black", markersize=8, label="CDS genomic span"),
    ]
    ax_len.legend(handles=legend_handles, fontsize=8, loc="upper right")

    # ══════════════════════════════════════════════════════════════════════════
    # Panel B: ORF Full-Span Rate (%)
    # ──────────────────────────────────────────────────────────────────────────
    # read가 CDS 전체 범위를 완전히 커버하는 비율입니다.
    # 비율이 높을수록 haplotype 정보가 풍부합니다.
    # 70% 이상이 권장되는 일반적인 기준입니다(nomadic2 기준).
    # 리드가 없는 유전자에 대해서는 NaN을 사용하고 그래프에서 제외합니다.
    # ══════════════════════════════════════════════════════════════════════════
    span_rates = []
    for gene in genes:
        sub = df[df["gene"] == gene]
        if sub.empty:
            span_rates.append(float("nan"))
        else:
            span_rates.append(100.0 * sub["spans_orf"].sum() / len(sub))

    # NaN 위치는 막대 그리지 않기 위해 nan-safe 처리
    bottoms = [0] * n_genes
    for i, (gene, rate) in enumerate(zip(genes, span_rates)):
        if not np.isnan(rate):
            bar = ax_span.bar(
                i, rate,
                color=colors[i], edgecolor="black", linewidth=0.5,
                zorder=2
            )
            # 막대 상단에 % 수치 표시
            ax_span.text(
                i, rate + 1.5, f"{rate:.0f}%",
                ha="center", va="bottom", fontsize=7
            )

    # 70% 권장 기준선 (빨간 점선)
    ax_span.axhline(70, color="red", lw=1.2, ls="--",
                    label="70% guideline", zorder=3)

    ax_span.set_xticks(range(n_genes))
    ax_span.set_xticklabels(genes, rotation=40, ha="right")
    ax_span.set_ylabel("Full-span reads (%)")
    ax_span.set_ylim(0, 115)
    ax_span.set_title("B. ORF Full-Span Rate", loc="left", fontweight="bold")
    ax_span.legend(fontsize=8)

    # ══════════════════════════════════════════════════════════════════════════
    # Panel C: AT Content (%) Boxplot
    # ──────────────────────────────────────────────────────────────────────────
    # 유전자별 AT 함량 분포를 boxplot으로 표시합니다.
    # Plasmodium falciparum 게놈의 전형적인 AT 함량(70-85%) 범위를
    # 연노란색 음영으로 표시합니다.
    # AT가 높을수록 homopolymer 오류와 GC-bias 관련 depth 이상이 생길 수 있습니다.
    # ══════════════════════════════════════════════════════════════════════════
    at_data = [df[df["gene"] == g]["at_pct"].dropna().values for g in genes]

    bp = ax_at.boxplot(
        at_data,
        positions=range(n_genes),
        widths=0.5,
        patch_artist=True,       # 박스 내부 색상 채우기 위해 필요
        flierprops={             # 이상치 점 스타일
            "marker": "o", "markersize": 3, "alpha": 0.4,
            "markerfacecolor": "gray", "markeredgecolor": "none"
        },
        medianprops={"color": "black", "lw": 1.5},  # 중앙선 강조
        whiskerprops={"lw": 0.8},
        capprops={"lw": 0.8},
    )
    # 박스 색상 유전자별로 지정
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.65)

    # Pf 전형 AT 범위 음영 (70–85%)
    ax_at.axhspan(70, 85, color="lightyellow", alpha=0.6,
                  label="Typical Pf range (70–85%)", zorder=0)
    ax_at.axhline(81, color="orange", lw=0.8, ls=":", label="Pf genome avg (~81%)")

    ax_at.set_xticks(range(n_genes))
    ax_at.set_xticklabels(genes, rotation=40, ha="right")
    ax_at.set_ylabel("AT content (%)")
    ax_at.set_title("C. AT Content per Gene", loc="left", fontweight="bold")
    ax_at.set_ylim(0, 105)
    ax_at.legend(fontsize=8)

    # ══════════════════════════════════════════════════════════════════════════
    # Panel D: Homopolymer Run Density (Scatter, Median per Gene)
    # ──────────────────────────────────────────────────────────────────────────
    # 각 유전자에 대해 ≥3, ≥6, ≥9 bp homopolymer run의
    # 중앙값(median) 밀도를 scatter plot으로 표시합니다.
    # 마커 모양으로 세 임계값을 구분합니다 (원, 사각형, 삼각형).
    # 각 유전자의 색상은 Panel A-C와 동일한 팔레트를 사용합니다.
    #
    # 해석:
    #   - ≥3 bp 밀도가 높은 유전자: 전체적으로 AT-rich 서열 구조
    #   - ≥9 bp 밀도가 높은 유전자: 긴 poly-A/T 반복이 많아 인델 오류 위험 높음
    # ══════════════════════════════════════════════════════════════════════════
    hp_specs = [
        ("hp3_pct", "≥3 bp",  "o",  60, 0.9),
        ("hp6_pct", "≥6 bp",  "s",  60, 0.75),
        ("hp9_pct", "≥9 bp",  "^",  60, 0.6),
    ]

    for col, label, marker, ms, alpha in hp_specs:
        medians = []
        for gene in genes:
            sub = df[df["gene"] == gene]
            if sub.empty or sub[col].isna().all():
                medians.append(float("nan"))
            else:
                medians.append(sub[col].median())

        ax_hp.scatter(
            range(n_genes), medians,
            marker=marker, s=ms,
            alpha=alpha, zorder=3,
            label=label,
            c=[colors[i] for i in range(n_genes)],  # 유전자별 색상 유지
            edgecolors="black", linewidths=0.5
        )

    ax_hp.set_xticks(range(n_genes))
    ax_hp.set_xticklabels(genes, rotation=40, ha="right")
    ax_hp.set_ylabel("Homopolymer density (% of read bases, median)")
    ax_hp.set_title("D. Homopolymer Run Density per Gene", loc="left", fontweight="bold")
    ax_hp.set_ylim(bottom=0)
    ax_hp.legend(fontsize=9, title="min. length", title_fontsize=8)

    # ── 그림 저장 ────────────────────────────────────────────────────────────
    plt.tight_layout()
    plt.savefig(out_pdf, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Figure saved: {out_pdf}", file=sys.stderr)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 5. CLI 및 메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--bam",        required=True,
                   help="정렬된 BAM 파일 (.bai 인덱스 필요)")
    p.add_argument("--bed",        required=True,
                   help="amplicons.bed (chrom start end gene)")
    p.add_argument("--cds",        required=True,
                   help="cds_coords.tsv")
    p.add_argument("--out-prefix", required=True, dest="out_prefix",
                   help="출력 파일 공통 접두사 (예: results/qc/barcode01_read_qc)")
    p.add_argument("--min-mapq",   type=int, default=0, dest="min_mapq",
                   help="최소 MAPQ 필터 (기본: 0)")
    p.add_argument("--min-qlen",   type=int, default=500, dest="min_qlen",
                   help="최소 read 길이 bp (기본: 500)")
    return p.parse_args()


def main():
    args = parse_args()

    _check_samtools()   # samtools PATH 확인을 main 진입 시점에 먼저 수행
    print(f"[INFO] BAM       : {args.bam}",        file=sys.stderr)
    print(f"[INFO] BED       : {args.bed}",        file=sys.stderr)
    print(f"[INFO] CDS       : {args.cds}",        file=sys.stderr)
    print(f"[INFO] min_mapq  : {args.min_mapq}",   file=sys.stderr)
    print(f"[INFO] min_qlen  : {args.min_qlen} bp", file=sys.stderr)

    # ── 좌표 로드 ────────────────────────────────────────────────────────────
    amplicons = load_amplicons(args.bed)
    cds_spans = load_cds_spans(args.cds)
    print(f"[INFO] Amplicons : {list(amplicons.keys())}", file=sys.stderr)

    # BAM 파일명에서 barcode 이름 추출
    # 예: "barcode01.Pf3D7.sorted.bam" → "barcode01"
    barcode = os.path.basename(args.bam).split(".")[0]

    # ── Read 통계 수집 ───────────────────────────────────────────────────────
    print("[INFO] Collecting read statistics from BAM...", file=sys.stderr)
    df = collect_read_stats(
        args.bam, amplicons, cds_spans,
        min_mapq=args.min_mapq, min_qlen=args.min_qlen
    )

    # 데이터가 없을 경우 빈 파일 생성 (파이프라인 중단 방지)
    if df.empty:
        print("[WARN] 수집된 read 없음. 빈 결과 파일을 생성합니다.", file=sys.stderr)
        df = pd.DataFrame(columns=[
            "gene", "read_len", "spans_orf",
            "at_pct", "hp3_pct", "hp6_pct", "hp9_pct"
        ])

    # ── TSV 저장 ─────────────────────────────────────────────────────────────
    os.makedirs(os.path.dirname(args.out_prefix) or ".", exist_ok=True)
    tsv_path = f"{args.out_prefix}_stats.tsv"
    df.to_csv(tsv_path, sep="\t", index=False)
    print(f"[INFO] Stats TSV : {tsv_path} ({len(df)} reads)", file=sys.stderr)

    # ── 콘솔 요약 출력 ───────────────────────────────────────────────────────
    if not df.empty:
        summary = (
            df.groupby("gene")
              .agg(
                  n_reads    = ("read_len",  "count"),
                  median_len = ("read_len",  "median"),
                  span_pct   = ("spans_orf", lambda x: round(100 * x.sum() / len(x), 1)),
                  median_at  = ("at_pct",    lambda x: round(x.median(), 1)),
                  median_hp3 = ("hp3_pct",   lambda x: round(x.median(), 1)),
              )
              .reset_index()
        )
        print("\n[SUMMARY]")
        print(summary.to_string(index=False))

    # ── 그림 생성 ────────────────────────────────────────────────────────────
    pdf_path = f"{args.out_prefix}.pdf"
    plot_read_qc(df, amplicons, cds_spans, pdf_path, barcode_label=barcode)


if __name__ == "__main__":
    main()
