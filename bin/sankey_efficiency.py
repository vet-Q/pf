#!/usr/bin/env python3
"""
sankey_efficiency.py
====================
배치(batch) 단위 시퀀싱 효율 Sankey 다이어그램 생성기.

nomadic2(Manser et al.) 논문 노트북 fig2_sequencing-efficiency-sankey의
분석 아이디어를 pfNextflow 파이프라인에 맞게 구현한 스크립트입니다.

── 시각화 흐름 ──────────────────────────────────────────────────────────────
  Total aligned reads
    ├── Unmapped               (Pf 게놈에 맵핑 실패 — 인간 DNA, 품질 불량 등)
    └── Mapped (Pf)
           ├── Off-target      (앰플리콘 영역 밖에 맵핑)
           └── On-target
                  ├── CRT
                  ├── DHFR
                  ├── DHPS
                  ├── K13
                  └── ... (BED 파일에 정의된 모든 유전자)

── 출처 대비 pfNextflow 차이점 ─────────────────────────────────────────────
  nomadic2는 인간 게놈(hs)에도 동시 정렬하여 Human DNA 오염 비율을 직접
  계산합니다. pfNextflow는 Pf 게놈에만 정렬하므로, "Unmapped" 비율이
  인간 DNA + 기타 미분류 reads의 합산으로 간주할 수 있습니다.
  임상 혈액 샘플에서 기생충혈증(parasitemia)이 낮으면 Unmapped 비율이
  높아지는 경향이 있습니다.

── 출력 ────────────────────────────────────────────────────────────────────
  sequencing_efficiency_sankey.html   Plotly 인터랙티브 HTML (마우스 오버 가능)
  sequencing_efficiency_sankey.pdf    정적 PDF 그림 (kaleido 필요)
  sequencing_efficiency_stats.tsv     barcode별 단계별 read 수 TSV

── 의존성 ──────────────────────────────────────────────────────────────────
  samtools     PATH에 설치되어 있어야 함
  pandas       pip install pandas
  plotly       pip install plotly
  kaleido      pip install kaleido   (PDF 저장용, 옵션)
               없으면 HTML만 저장됩니다.

── 사용 예 (단독 실행) ─────────────────────────────────────────────────────
  # BAM 목록 파일 방식
  python3 bin/sankey_efficiency.py \\
      --bam-list bam_list.txt \\
      --bed      resources/amplicons.bed \\
      --out-dir  results/qc/

  # BAM 디렉토리 방식 (*.sorted.bam 자동 검색)
  python3 bin/sankey_efficiency.py \\
      --bam-dir  results/alignment/minimap2/ \\
      --bed      resources/amplicons.bed \\
      --out-dir  results/qc/
"""

import argparse
import glob
import os
import subprocess
import sys

import pandas as pd

# plotly는 선택적 의존성 — 없으면 오류 메시지를 명확히 출력
try:
    import plotly.graph_objects as go
    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False

# kaleido는 PDF 저장용 선택적 의존성
try:
    import kaleido  # noqa: F401 (import 성공 여부만 확인)
    HAS_KALEIDO = True
except ImportError:
    HAS_KALEIDO = False


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 1. 좌표 파일 로드
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def load_amplicons(bed_path):
    """
    amplicons.bed를 읽어 [(gene, chrom, start, end), ...] 리스트를 반환합니다.
    gene별로 samtools view -c region 명령에 사용할 좌표를 가져옵니다.

    반환: [ ("CRT", "Pf3D7_07_v3", 402384, 406341), ... ]
    """
    genes = []
    with open(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            chrom = parts[0]
            start = int(parts[1])   # 0-based
            end   = int(parts[2])   # 0-based exclusive
            gene  = parts[3]
            genes.append((gene, chrom, start, end))
    return genes


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 2. samtools를 이용한 read count 수집
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def run_cmd(cmd):
    """
    shell 명령을 실행하고 stdout 문자열을 반환합니다.
    실패 시 오류 메시지를 출력하고 빈 문자열 반환 (파이프라인 중단 방지).
    """
    try:
        result = subprocess.run(
            cmd, shell=True,
            capture_output=True, text=True, check=True
        )
        return result.stdout.strip()
    except subprocess.CalledProcessError as e:
        print(f"[WARN] 명령 실패: {cmd}\n       {e.stderr.strip()}", file=sys.stderr)
        return ""


def get_flagstat_counts(bam_path):
    """
    samtools flagstat로 BAM 파일의 기본 통계를 파싱합니다.

    반환: dict { "total": int, "mapped": int }

    flagstat 출력 예:
        10000 + 0 in total (QC-passed reads + QC-failed reads)
        ...
        9500 + 0 mapped (95.00% : N/A)
        ...

    주의:
        flagstat의 "total"에는 supplementary alignment가 포함됩니다.
        Nanopore long-read 데이터에서는 split alignment가 있어
        실제 read 수보다 클 수 있습니다.
        "-F 2048"(supplementary 제외) 플래그를 사용하면 더 정확합니다.
    """
    # supplementary (-F 2048) 및 secondary (-F 256) 제외한 기본 reads만 카운트
    total_raw = run_cmd(f"samtools view -c -F 2304 {bam_path}")
    mapped_raw = run_cmd(f"samtools view -c -F 2308 {bam_path}")
    # -F 2308 = -F (2048 supplementary + 256 secondary + 4 unmapped)

    try:
        total  = int(total_raw)
        mapped = int(mapped_raw)
    except ValueError:
        print(f"[WARN] {bam_path}: flagstat 파싱 실패, 0으로 설정", file=sys.stderr)
        total = mapped = 0

    unmapped = total - mapped
    return {"total": total, "mapped": mapped, "unmapped": unmapped}


def get_ontarget_counts(bam_path, bed_path):
    """
    samtools view -c -L <BED>로 앰플리콘 영역 내 reads 수를 셉니다.

    "-L bed" 옵션: BED 파일의 모든 영역과 겹치는 reads를 카운트합니다.
    supplementary(-F 2048), secondary(-F 256) 제외합니다.

    반환: int (on-target read 수)
    """
    count_raw = run_cmd(
        f"samtools view -c -F 2304 -L {bed_path} {bam_path}"
    )
    try:
        return int(count_raw)
    except ValueError:
        return 0


def get_gene_counts(bam_path, amplicons):
    """
    각 앰플리콘/유전자 region의 read 수를 samtools view -c로 셉니다.

    samtools 좌표: "chrom:start+1-end" (1-based, 포함형)
    BED 좌표(0-based half-open) → samtools 좌표(1-based closed) 변환:
        samtools_start = bed_start + 1
        samtools_end   = bed_end   (= 0-based exclusive와 동일, 1-based inclusive)

    반환: { "CRT": 1234, "DHFR": 567, ... }
    """
    counts = {}
    for gene, chrom, bed_start, bed_end in amplicons:
        # samtools region 문자열 생성 (1-based)
        region = f"{chrom}:{bed_start + 1}-{bed_end}"
        count_raw = run_cmd(
            f"samtools view -c -F 2304 {bam_path} {region}"
        )
        try:
            counts[gene] = int(count_raw)
        except ValueError:
            counts[gene] = 0
    return counts


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 3. 전체 배치 집계
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def collect_batch_stats(bam_paths, bed_path, amplicons):
    """
    모든 BAM 파일에서 단계별 read count를 수집합니다.

    처리 순서:
      1. flagstat → total, mapped, unmapped
      2. samtools view -L bed → on-target 수
         off-target = mapped - on-target
      3. 각 유전자 region → gene별 수

    반환:
        all_rows  : 바코드별 상세 DataFrame
        batch_sum : 배치 전체 합산 Series
    """
    gene_names = [g[0] for g in amplicons]
    rows = []

    for bam_path in bam_paths:
        barcode = os.path.basename(bam_path).split(".")[0]
        print(f"[INFO] Processing {barcode} ...", file=sys.stderr)

        # ── flagstat ─────────────────────────────────────────────────────────
        fc = get_flagstat_counts(bam_path)

        # ── on-target / off-target ────────────────────────────────────────────
        on_target  = get_ontarget_counts(bam_path, bed_path)
        # off_target: mapped 중 on-target이 아닌 수
        # (supplementary를 제외했으므로 단순 뺄셈으로 계산)
        off_target = max(0, fc["mapped"] - on_target)

        # ── 유전자별 count ────────────────────────────────────────────────────
        gene_counts = get_gene_counts(bam_path, amplicons)

        row = {
            "barcode"   : barcode,
            "total"     : fc["total"],
            "mapped"    : fc["mapped"],
            "unmapped"  : fc["unmapped"],
            "on_target" : on_target,
            "off_target": off_target,
        }
        # 유전자별 컬럼 추가
        for gene in gene_names:
            row[gene] = gene_counts.get(gene, 0)

        rows.append(row)

    df        = pd.DataFrame(rows)
    batch_sum = df.drop(columns="barcode").sum()
    return df, batch_sum


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 4. Plotly Sankey 다이어그램 생성
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def make_sankey(batch_sum, gene_names, title="Sequencing Efficiency"):
    """
    Plotly Sankey 다이어그램 Figure를 생성합니다.

    노드 구성:
        0: Total aligned reads
        1: Unmapped
        2: Mapped (Pf)
        3: Off-target
        4: On-target
        5 ~ 5+n: 유전자별 (CRT, DHFR, ...)

    링크(edge) 구성:
        Total  → Unmapped    (batch_sum["unmapped"])
        Total  → Mapped      (batch_sum["mapped"])
        Mapped → Off-target  (batch_sum["off_target"])
        Mapped → On-target   (batch_sum["on_target"])
        On-target → Gene_i   (batch_sum[gene_i])

    각 링크에는 read 수와 비율(%)이 호버 텍스트로 표시됩니다.

    매개변수:
        batch_sum  : 배치 합산 pd.Series
        gene_names : 유전자 이름 목록 (BED 파일 순서)
        title      : 그림 제목

    반환:
        plotly.graph_objects.Figure
    """
    total      = max(1, int(batch_sum.get("total",      0)))
    unmapped   = int(batch_sum.get("unmapped",   0))
    mapped     = int(batch_sum.get("mapped",     0))
    off_target = int(batch_sum.get("off_target", 0))
    on_target  = int(batch_sum.get("on_target",  0))

    # ── 퍼센트 계산 헬퍼 ────────────────────────────────────────────────────
    def pct(val, denom):
        """val / denom * 100, denom=0이면 0 반환"""
        return round(100.0 * val / denom, 1) if denom > 0 else 0.0

    # ── 노드 정의 ────────────────────────────────────────────────────────────
    # 인덱스 순서가 links의 source/target과 대응되어야 합니다.
    node_labels = (
        [
            f"Total\n{total:,}",
            f"Unmapped\n{unmapped:,}\n({pct(unmapped, total)}%)",
            f"Mapped (Pf)\n{mapped:,}\n({pct(mapped, total)}%)",
            f"Off-target\n{off_target:,}\n({pct(off_target, mapped)}% of mapped)",
            f"On-target\n{on_target:,}\n({pct(on_target, mapped)}% of mapped)",
        ]
        + [
            # 유전자별 노드: on-target 내 비율 표시
            f"{gene}\n{int(batch_sum.get(gene, 0)):,}\n"
            f"({pct(int(batch_sum.get(gene, 0)), on_target)}% of on-target)"
            for gene in gene_names
        ]
    )

    # 노드 색상 팔레트 (Plotly 내장 색상 이름 사용)
    node_colors = (
        ["#2196F3",          # Total — 파란색
         "#9E9E9E",          # Unmapped — 회색
         "#4CAF50",          # Mapped (Pf) — 초록색
         "#FF9800",          # Off-target — 주황색
         "#3F51B5"]          # On-target — 남색
        + [                  # 유전자별 — teal 계열 순환
            f"hsl({200 + i * 25}, 60%, 50%)"
            for i in range(len(gene_names))
        ]
    )

    # ── 링크 정의 (source, target, value) ────────────────────────────────────
    # 인덱스 매핑:
    #   0=Total, 1=Unmapped, 2=Mapped, 3=Off-target, 4=On-target,
    #   5+i=gene_i

    link_sources  = []
    link_targets  = []
    link_values   = []
    link_labels   = []   # 호버 텍스트
    link_colors   = []

    def add_link(src, tgt, val, label, color="rgba(150,150,150,0.35)"):
        """링크 추가 헬퍼 — val이 0이면 추가하지 않음 (Sankey 시각적 오류 방지)"""
        if val > 0:
            link_sources.append(src)
            link_targets.append(tgt)
            link_values.append(val)
            link_labels.append(label)
            link_colors.append(color)

    # Total → Unmapped
    add_link(0, 1, unmapped,
             f"Unmapped: {unmapped:,} ({pct(unmapped, total)}%)",
             "rgba(158,158,158,0.4)")

    # Total → Mapped (Pf)
    add_link(0, 2, mapped,
             f"Mapped to Pf: {mapped:,} ({pct(mapped, total)}%)",
             "rgba(76,175,80,0.4)")

    # Mapped → Off-target
    add_link(2, 3, off_target,
             f"Off-target: {off_target:,} ({pct(off_target, mapped)}%)",
             "rgba(255,152,0,0.4)")

    # Mapped → On-target
    add_link(2, 4, on_target,
             f"On-target: {on_target:,} ({pct(on_target, mapped)}%)",
             "rgba(63,81,181,0.4)")

    # On-target → 각 유전자
    for i, gene in enumerate(gene_names):
        gene_count = int(batch_sum.get(gene, 0))
        add_link(
            4,      # On-target 노드 인덱스
            5 + i,  # 해당 유전자 노드 인덱스
            gene_count,
            f"{gene}: {gene_count:,} ({pct(gene_count, on_target)}% of on-target)",
            f"hsla({200 + i * 25}, 60%, 50%, 0.45)"
        )

    # ── Sankey Figure 생성 ───────────────────────────────────────────────────
    fig = go.Figure(go.Sankey(
        arrangement="snap",   # 노드를 자동 배치 후 열 기준으로 snap
        node=dict(
            pad          = 20,          # 노드 간 수직 간격(px)
            thickness    = 25,          # 노드 너비(px)
            line         = dict(color="black", width=0.5),
            label        = node_labels,
            color        = node_colors,
            hovertemplate= "%{label}<extra></extra>",
        ),
        link=dict(
            source         = link_sources,
            target         = link_targets,
            value          = link_values,
            label          = link_labels,
            color          = link_colors,
            hovertemplate  = "%{label}<extra></extra>",
        ),
    ))

    fig.update_layout(
        title=dict(
            text  = title,
            font  = dict(size=16),
            x     = 0.5,
            xanchor = "center",
        ),
        font      = dict(size=11),
        height    = 600,
        width     = 1000,
        margin    = dict(l=20, r=20, t=60, b=20),
    )
    return fig


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 5. barcode별 비교 히트맵 (보조 출력)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def make_barcode_heatmap(df, gene_names):
    """
    barcode별 주요 지표를 Plotly 히트맵으로 시각화합니다.

    표시 지표 (열):
        - mapped_pct   (%) : Pf 게놈 맵핑 성공 비율
        - ontarget_pct (%) : mapped 중 on-target 비율
        - gene별          : on-target 중 각 유전자 비율

    히트맵에서 발견할 수 있는 패턴:
        - mapped_pct 낮은 barcode: 기생충혈증 낮거나 DNA 품질 불량
        - ontarget_pct 낮은 barcode: 비특이적 증폭 문제
        - 특정 유전자 열이 낮은 barcode: 해당 amplicon의 선택적 실패

    반환: plotly.graph_objects.Figure
    """
    # 백분율 계산
    metrics = df.copy()
    metrics["mapped_pct"]   = (100 * metrics["mapped"]    / metrics["total"].clip(1)).round(1)
    metrics["ontarget_pct"] = (100 * metrics["on_target"] / metrics["mapped"].clip(1)).round(1)
    for gene in gene_names:
        metrics[f"{gene}_pct"] = (
            100 * metrics[gene] / metrics["on_target"].clip(1)
        ).round(1)

    # 히트맵용 행렬 구성
    col_keys    = ["mapped_pct", "ontarget_pct"] + [f"{g}_pct" for g in gene_names]
    col_labels  = ["Mapped (%)", "On-target (%)"] + [f"{g} (%)" for g in gene_names]
    z_matrix    = metrics[col_keys].values.T.tolist()   # (n_cols, n_barcodes)
    barcodes    = metrics["barcode"].tolist()

    fig = go.Figure(go.Heatmap(
        z             = z_matrix,
        x             = barcodes,
        y             = col_labels,
        colorscale    = "RdYlGn",   # 낮음=빨강, 높음=초록으로 직관적 시각화
        zmin          = 0,
        zmax          = 100,
        text          = [[f"{v:.0f}%" for v in row] for row in z_matrix],
        texttemplate  = "%{text}",
        colorbar      = dict(title="Percentage"),
        hovertemplate = "Barcode: %{x}<br>Metric: %{y}<br>Value: %{z:.1f}%<extra></extra>",
    ))

    fig.update_layout(
        title=dict(
            text    = "Per-Barcode Sequencing Efficiency",
            font    = dict(size=14),
            x       = 0.5,
            xanchor = "center",
        ),
        xaxis=dict(tickangle=-45, title="Barcode"),
        yaxis=dict(title="Metric"),
        height = 400 + len(col_labels) * 20,
        width  = 200 + len(barcodes) * 60,
        margin = dict(l=120, r=20, t=60, b=100),
    )
    return fig


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 6. CLI 및 메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # BAM 파일 지정 방법: 목록 파일 또는 디렉토리 중 하나 필수
    grp = p.add_mutually_exclusive_group(required=True)
    grp.add_argument(
        "--bam-list", dest="bam_list",
        help="BAM 파일 경로를 한 줄씩 기록한 텍스트 파일"
    )
    grp.add_argument(
        "--bam-dir", dest="bam_dir",
        help="BAM 파일이 있는 디렉토리 (*.sorted.bam 자동 검색)"
    )

    p.add_argument("--bed",     required=True, help="amplicons.bed")
    p.add_argument("--out-dir", required=True, dest="out_dir",
                   help="출력 디렉토리 경로")
    p.add_argument("--title",   default="Sequencing Efficiency — pfNextflow",
                   help="그림 제목 (기본: 'Sequencing Efficiency — pfNextflow')")
    return p.parse_args()


def main():
    args = parse_args()

    # plotly 없으면 즉시 오류 종료
    if not HAS_PLOTLY:
        sys.exit(
            "[ERROR] plotly 패키지가 없습니다.\n"
            "        pip install plotly  또는\n"
            "        conda install -c conda-forge plotly  으로 설치하세요."
        )

    # ── BAM 파일 목록 수집 ────────────────────────────────────────────────────
    if args.bam_list:
        with open(args.bam_list) as fh:
            bam_paths = [line.strip() for line in fh if line.strip()]
    else:
        # 디렉토리에서 *.sorted.bam 파일 자동 검색
        bam_paths = sorted(glob.glob(os.path.join(args.bam_dir, "*.sorted.bam")))
        if not bam_paths:
            # *.bam 패턴도 시도
            bam_paths = sorted(glob.glob(os.path.join(args.bam_dir, "*.bam")))

    if not bam_paths:
        sys.exit("[ERROR] BAM 파일을 찾을 수 없습니다.")

    # BAI 인덱스 파일 존재 여부 확인 (없으면 경고)
    for bp in bam_paths:
        if not os.path.exists(bp + ".bai") and not os.path.exists(bp[:-4] + ".bai"):
            print(f"[WARN] BAI 인덱스 없음: {bp} — samtools index {bp} 실행 필요",
                  file=sys.stderr)

    print(f"[INFO] {len(bam_paths)}개 BAM 파일 발견", file=sys.stderr)

    # ── amplicons.bed 파싱 ────────────────────────────────────────────────────
    amplicons  = load_amplicons(args.bed)
    gene_names = [g[0] for g in amplicons]
    print(f"[INFO] Amplicons : {gene_names}", file=sys.stderr)

    # ── 출력 디렉토리 생성 ────────────────────────────────────────────────────
    os.makedirs(args.out_dir, exist_ok=True)

    # ── 배치 통계 수집 ────────────────────────────────────────────────────────
    print("[INFO] Collecting read statistics (this may take a few minutes)...",
          file=sys.stderr)
    df, batch_sum = collect_batch_stats(bam_paths, args.bed, amplicons)

    # ── TSV 저장 ─────────────────────────────────────────────────────────────
    tsv_path = os.path.join(args.out_dir, "sequencing_efficiency_stats.tsv")
    df.to_csv(tsv_path, sep="\t", index=False)
    print(f"[INFO] Stats TSV: {tsv_path}", file=sys.stderr)

    # ── 배치 합산 콘솔 요약 ──────────────────────────────────────────────────
    total   = max(1, int(batch_sum.get("total", 0)))
    mapped  = int(batch_sum.get("mapped", 0))
    on_tgt  = int(batch_sum.get("on_target", 0))
    print("\n[BATCH SUMMARY]")
    print(f"  Total    : {total:>10,}")
    print(f"  Mapped   : {mapped:>10,}  ({100*mapped/total:.1f}%)")
    print(f"  On-target: {on_tgt:>10,}  ({100*on_tgt/mapped:.1f}% of mapped)" if mapped else "")
    for gene in gene_names:
        g_cnt = int(batch_sum.get(gene, 0))
        print(f"  {gene:<10}: {g_cnt:>10,}  ({100*g_cnt/on_tgt:.1f}% of on-target)" if on_tgt else "")

    # ── Sankey 그림 생성 ─────────────────────────────────────────────────────
    sankey_fig   = make_sankey(batch_sum, gene_names, title=args.title)
    heatmap_fig  = make_barcode_heatmap(df, gene_names)

    # HTML 저장 (항상 가능 — plotly만 있으면 됨)
    html_path = os.path.join(args.out_dir, "sequencing_efficiency_sankey.html")
    sankey_fig.write_html(html_path)
    print(f"[INFO] Sankey HTML: {html_path}", file=sys.stderr)

    heatmap_html = os.path.join(args.out_dir, "barcode_efficiency_heatmap.html")
    heatmap_fig.write_html(heatmap_html)
    print(f"[INFO] Heatmap HTML: {heatmap_html}", file=sys.stderr)

    # PDF 저장 (kaleido 패키지가 있을 때만)
    if HAS_KALEIDO:
        pdf_path = os.path.join(args.out_dir, "sequencing_efficiency_sankey.pdf")
        sankey_fig.write_image(pdf_path, format="pdf", width=1000, height=600)
        print(f"[INFO] Sankey PDF: {pdf_path}", file=sys.stderr)

        heatmap_pdf = os.path.join(args.out_dir, "barcode_efficiency_heatmap.pdf")
        heatmap_fig.write_image(heatmap_pdf, format="pdf",
                                width=200 + len(bam_paths) * 60,
                                height=400 + len(gene_names) * 25)
        print(f"[INFO] Heatmap PDF: {heatmap_pdf}", file=sys.stderr)
    else:
        print(
            "[INFO] kaleido 미설치 → PDF 생략. HTML 파일을 브라우저에서 확인하세요.\n"
            "       pip install kaleido  으로 설치하면 PDF도 생성됩니다.",
            file=sys.stderr
        )

    print("[INFO] Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
