#!/usr/bin/env python3
"""
splice_and_translate.py

Two modes:

  Consensus mode (default, --fasta):
    Reads a per-sample consensus FASTA produced by bcftools consensus.
    The FASTA header carries the samtools region (e.g. Pf3D7_07_v3:402385-406341),
    which is used to calculate the offset between genomic coordinates and
    sequence positions.

  Reference mode (--ref-fasta):
    Reads from the reference genome FASTA directly via samtools faidx.
    Used once per gene to generate the 3D7 reference sequence for alignment.

In both modes the script:
  1. Loads CDS exon boundaries from cds_coords.tsv
  2. Extracts each exon's sequence from the input
  3. Joins exons in transcript order (ascending for +, descending for - strand)
  4. Reverse-complements the joined sequence for minus-strand genes
  5. Translates to protein (standard codon table, Ns → X)
  6. Writes one-record FASTA files for nucleotide and amino acid output

Usage:
  # Consensus mode
  splice_and_translate.py \\
      --fasta   barcode01.CRT.bcftools.consensus.fasta \\
      --cds     cds_coords.tsv \\
      --gene    CRT \\
      --sample  barcode01 \\
      --out-nt  barcode01.CRT.nt.fasta \\
      --out-aa  barcode01.CRT.aa.fasta

  # Reference mode
  splice_and_translate.py \\
      --ref-fasta PlasmoDB-67_Pfalciparum3D7_Genome.fasta \\
      --cds       cds_coords.tsv \\
      --gene      CRT \\
      --sample    3D7 \\
      --out-nt    3D7.CRT.nt.fasta \\
      --out-aa    3D7.CRT.aa.fasta
"""

import argparse
import gzip
import os
import subprocess
import sys

# ── Standard genetic code ────────────────────────────────────────────────────
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

_RC = str.maketrans(
    'ACGTacgtNnRrYyKkMmSsWwBbDdHhVv',
    'TGCAtgcaNnYyRrMmKkSsWwVvHhDdBb'
)

# bcftools consensus 가 heterozygous site 에 삽입하는 IUPAC ambiguity 코드와
# 그것이 나타낼 수 있는 표준 염기 집합의 매핑.
# VCF 의 AF 를 이용해 어느 쪽 염기인지 결정할 때 사용한다.
IUPAC_BASES: dict = {
    'R': frozenset('AG'),  'Y': frozenset('CT'),  'S': frozenset('GC'),
    'W': frozenset('AT'),  'K': frozenset('GT'),  'M': frozenset('AC'),
    'B': frozenset('CGT'), 'D': frozenset('AGT'),
    'H': frozenset('ACT'), 'V': frozenset('ACG'),
    'N': frozenset('ACGT'),
}


# ── Sequence utilities ────────────────────────────────────────────────────────

def rev_comp(seq: str) -> str:
    return seq.translate(_RC)[::-1]


def translate(seq: str) -> str:
    """Translate nucleotide sequence to protein. Ambiguous codons → X."""
    aa = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3].upper()
        if len(codon) < 3:
            break
        if any(c not in 'ACGT' for c in codon):
            aa.append('X')
        else:
            aa.append(CODON_TABLE.get(codon, 'X'))
    return ''.join(aa)


# ── I/O helpers ───────────────────────────────────────────────────────────────

def read_fasta_first(path: str):
    """Read the first record from a FASTA file. Returns (header, seq_uppercase)."""
    header, parts = '', []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if header:
                    break
                header = line[1:]
            else:
                parts.append(line.upper())
    return header, ''.join(parts)


def faidx_extract(ref_fasta: str, chrom: str, start0: int, end0: int):
    """
    Extract a sequence from a FASTA using samtools faidx.
    start0  : 0-based inclusive
    end0    : 0-based exclusive
    Returns : (region_string, seq_uppercase)
    """
    region = f"{chrom}:{start0 + 1}-{end0}"   # samtools uses 1-based, inclusive
    result = subprocess.run(
        ['samtools', 'faidx', ref_fasta, region],
        capture_output=True, text=True, check=True,
    )
    lines = result.stdout.strip().split('\n')
    if not lines or not lines[0].startswith('>'):
        raise RuntimeError(f"samtools faidx returned unexpected output for {region}")
    seq = ''.join(lines[1:]).upper()
    return region, seq


def get_bed_start(header: str) -> int:
    """
    Parse the 0-based genome start position from a samtools FASTA header.
    Example header: 'Pf3D7_07_v3:402385-406341'
    samtools outputs 1-based coords → convert to 0-based by subtracting 1.
    """
    try:
        region  = header.split()[0]          # strip any trailing description
        start1  = int(region.split(':')[1].split('-')[0])
        return start1 - 1
    except Exception:
        return 0


# ── CDS loader ────────────────────────────────────────────────────────────────

def load_cds(tsv: str, gene: str) -> list:
    """
    Return a list of exon dicts for *gene* from cds_coords.tsv.
    Each dict: {chrom, start (0-based), end (0-based excl), strand, rank}
    Sorted by genomic start ascending.
    """
    exons = []
    with open(tsv) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 6 or parts[0] != gene:
                continue
            exons.append({
                'chrom':  parts[1],
                'start':  int(parts[2]),
                'end':    int(parts[3]),
                'strand': parts[4],
                'rank':   int(parts[5]),
            })
    return sorted(exons, key=lambda x: x['start'])


# ── Splicing ──────────────────────────────────────────────────────────────────

def splice_exons(seq: str, seq_start0: int, exons: list, strand: str) -> str:
    """
    Extract exon sub-sequences from *seq* and join them in transcript order.

    For + strand : join exons in ascending genomic order.
    For - strand : join exons in DESCENDING genomic order (5'→3' of transcript),
                   then reverse-complement the concatenated result.

    Parameters
    ----------
    seq        : nucleotide string covering the region (0-based position 0
                 corresponds to genomic coordinate seq_start0)
    seq_start0 : 0-based genomic coordinate of position 0 in *seq*
    exons      : list of exon dicts (already sorted by start ascending)
    strand     : '+' or '-'
    """
    ordered = exons if strand == '+' else list(reversed(exons))

    parts = []
    for ex in ordered:
        o_start = ex['start'] - seq_start0
        o_end   = ex['end']   - seq_start0
        if o_start < 0 or o_end > len(seq):
            print(
                f"  WARNING: exon rank={ex['rank']} offset [{o_start}:{o_end}] "
                f"out of bounds (seq len={len(seq)})",
                file=sys.stderr,
            )
            o_start = max(0, o_start)
            o_end   = min(len(seq), o_end)
        parts.append(seq[o_start:o_end])

    joined = ''.join(parts)
    return rev_comp(joined) if strand == '-' else joined


# ── VCF AF parser ─────────────────────────────────────────────────────────────

def parse_vcf_af(vcf_path: str) -> dict:
    """
    bcftools VCF(gz)에서 SNP 위치별 allele frequency 를 읽는다.
    FORMAT 에 AD(allele depth) 가 있어야 하며,
    AF = alt_AD / (ref_AD + alt_AD) 로 계산한다.

    Returns
    -------
    dict : { genomic_pos_0based : {
               'chrom': str, 'ref': str, 'alt': str,
               'ref_ad': int, 'alt_ad': int, 'af': float } }
    """
    import gzip as _gzip
    lookup: dict = {}
    opener = _gzip.open if vcf_path.endswith('.gz') else open
    with opener(vcf_path, 'rt') as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) < 10:
                continue
            chrom, pos1, ref, alt_field = parts[0], int(parts[1]), parts[3], parts[4]
            alt = alt_field.split(',')[0]          # 첫 번째 ALT allele 만
            if len(ref) != 1 or len(alt) != 1:
                continue                            # indel 제외
            fmt = dict(zip(parts[8].split(':'), parts[9].rstrip().split(':')))
            ad_str = fmt.get('AD', '.')
            if '.' in ad_str:
                continue
            ads = ad_str.split(',')
            if len(ads) < 2:
                continue
            try:
                ref_ad, alt_ad = int(ads[0]), int(ads[1])
            except ValueError:
                continue
            total = ref_ad + alt_ad
            if total == 0:
                continue
            lookup[pos1 - 1] = {               # 0-based 변환
                'chrom': chrom, 'ref': ref, 'alt': alt,
                'ref_ad': ref_ad, 'alt_ad': alt_ad, 'af': alt_ad / total,
            }
    return lookup


# ── CDS position map ──────────────────────────────────────────────────────────

def build_cds_pos_map(exons: list, strand: str) -> dict:
    """
    게놈 위치(0-based) → CDS 내 위치(0-based) 매핑을 반환한다.

    + strand : exon 을 오름차순 처리, exon 내부도 5'→3' 방향.
    - strand : exon 을 내림차순 처리(높은 게놈 좌표가 transcript 5'),
               exon 내부도 높은 쪽에서 낮은 쪽으로 순회.
    """
    pos_map: dict = {}
    cds_i = 0
    if strand == '+':
        for ex in exons:
            for gpos in range(ex['start'], ex['end']):
                pos_map[gpos] = cds_i
                cds_i += 1
    else:
        for ex in reversed(exons):
            for gpos in range(ex['end'] - 1, ex['start'] - 1, -1):
                pos_map[gpos] = cds_i
                cds_i += 1
    return pos_map


# ── IUPAC 해소 ── VCF AF 기반 ─────────────────────────────────────────────────

def resolve_and_report(
    seq: str,
    seq_start0: int,
    vcf_lookup: dict,
    cds_pos_map: dict,
    gene: str,
    sample: str,
    threshold: float,
) -> tuple:
    """
    consensus FASTA 의 IUPAC 코드를 VCF AF 로 해소한다.

    결정 규칙
    ---------
    af >= threshold          → ALT 염기로 콜  ('alt_call')
    af <= 1 - threshold      → REF 염기로 콜  ('ref_call')
    그 외 (AF 가 중간 범위)    → IUPAC 유지   ('ambiguous_mixed') ← 혼합감염 의심

    Returns
    -------
    resolved_seq : str   — IUPAC 해소된 서열
    af_records   : list  — CDS 내 IUPAC 위치마다 한 행, TSV 출력에 사용
    """
    chars = list(seq)
    af_records: list = []

    for i, base in enumerate(chars):
        base_up = base.upper()
        if base_up not in IUPAC_BASES:
            continue                            # 일반 ACGT — 건너뜀

        gpos    = seq_start0 + i               # 0-based 게놈 좌표
        cds_pos = cds_pos_map.get(gpos)        # CDS 밖이면 None

        if gpos not in vcf_lookup:
            continue                            # VCF 레코드 없음 — 건너뜀

        info = vcf_lookup[gpos]
        ref, alt  = info['ref'], info['alt']
        ref_ad, alt_ad, af = info['ref_ad'], info['alt_ad'], info['af']

        if af >= threshold:
            resolved = alt
            note     = 'alt_call'
        elif af <= 1.0 - threshold:
            resolved = ref
            note     = 'ref_call'
        else:
            resolved = base_up                 # IUPAC 유지 (혼합감염 의심)
            note     = 'ambiguous_mixed'

        chars[i] = resolved

        if cds_pos is not None:
            codon_num   = cds_pos // 3 + 1    # 1-based 코돈 번호
            codon_frame = cds_pos % 3          # 0=1st, 1=2nd, 2=3rd position
            af_records.append({
                'sample':      sample,
                'gene':        gene,
                'chrom':       info['chrom'],
                'gen_pos':     gpos + 1,       # 1-based (보고용)
                'cds_pos':     cds_pos + 1,    # 1-based
                'codon':       codon_num,
                'frame':       codon_frame,
                'iupac':       base_up,
                'ref':         ref,
                'alt':         alt,
                'ref_depth':   ref_ad,
                'alt_depth':   alt_ad,
                'total_depth': ref_ad + alt_ad,
                'af':          round(af, 4),
                'resolved':    resolved,
                'note':        note,
            })

    return ''.join(chars), af_records


# ── AF 표 저장 ────────────────────────────────────────────────────────────────

def write_af_table(records: list, path: str) -> None:
    """
    IUPAC 해소 정보를 TSV 로 저장한다.
    records 가 비어있어도 헤더만 있는 파일을 생성한다 (Nextflow 패턴 매칭용).
    """
    cols = [
        'sample', 'gene', 'chrom', 'gen_pos', 'cds_pos', 'codon', 'frame',
        'iupac', 'ref', 'alt', 'ref_depth', 'alt_depth', 'total_depth',
        'af', 'resolved', 'note',
    ]
    with open(path, 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for rec in records:
            fh.write('\t'.join(str(rec[c]) for c in cols) + '\n')


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description='Splice CDS exons from a consensus or reference FASTA and translate'
    )
    ap.add_argument(
        '--fasta',
        help='Per-sample consensus FASTA (bcftools consensus output) [consensus mode]',
    )
    ap.add_argument(
        '--ref-fasta',
        help='Reference genome FASTA (accessed via samtools faidx) [reference mode]',
    )
    ap.add_argument('--cds',    required=True, help='CDS coordinates TSV (cds_coords.tsv)')
    ap.add_argument('--gene',   required=True, help='Gene name (e.g. CRT)')
    ap.add_argument('--sample', required=True, help='Sample label written to FASTA header')
    ap.add_argument('--out-nt', required=True, help='Output nucleotide FASTA path')
    ap.add_argument('--out-aa', required=True, help='Output amino-acid FASTA path')
    ap.add_argument(
        '--vcf',
        default=None,
        help='Per-gene bcftools VCF (.vcf.gz, with FORMAT/AD) for IUPAC resolution'
    )
    ap.add_argument(
        '--af-threshold', type=float, default=0.6,
        help='AF >= t → call ALT; AF <= 1-t → call REF; otherwise keep IUPAC (default 0.6)'
    )
    ap.add_argument(
        '--out-af', default=None,
        help='Output TSV path for IUPAC AF resolution table'
    )
    args = ap.parse_args()

    if not args.fasta and not args.ref_fasta:
        ap.error('Provide either --fasta (consensus mode) or --ref-fasta (reference mode)')

    # ── Load CDS definitions ──────────────────────────────────────────────────
    exons = load_cds(args.cds, args.gene)
    if not exons:
        print(f"WARNING: no CDS entries for gene '{args.gene}' in {args.cds}",
              file=sys.stderr)
        open(args.out_nt, 'w').close()
        open(args.out_aa, 'w').close()
        return

    strand = exons[0]['strand']

    # ── Obtain sequence and its genomic anchor position ───────────────────────
    if args.ref_fasta:
        # Reference mode: extract minimal spanning region from genome
        chrom   = exons[0]['chrom']
        span_s  = min(e['start'] for e in exons)
        span_e  = max(e['end']   for e in exons)
        _hdr, seq = faidx_extract(args.ref_fasta, chrom, span_s, span_e)
        seq_start0 = span_s
    else:
        # Consensus mode: read FASTA and parse header for offset
        header, seq = read_fasta_first(args.fasta)
        seq_start0  = get_bed_start(header)

    # ── IUPAC 해소 (consensus 모드 + VCF 제공 시) ──────────────────────────────
    af_records: list = []
    if args.fasta and args.vcf:
        if not os.path.exists(args.vcf):
            print(f"  WARNING: --vcf not found, skipping IUPAC resolution: {args.vcf}",
                  file=sys.stderr)
        else:
            vcf_lookup = parse_vcf_af(args.vcf)
            if vcf_lookup:
                cds_pos_map = build_cds_pos_map(exons, strand)
                seq, af_records = resolve_and_report(
                    seq, seq_start0, vcf_lookup, cds_pos_map,
                    args.gene, args.sample, args.af_threshold,
                )
                n_resolved = sum(1 for r in af_records if r['note'] != 'ambiguous_mixed')
                n_mixed    = len(af_records) - n_resolved
                if af_records:
                    print(
                        f"  {args.gene} [{args.sample}]: "
                        f"{n_resolved} IUPAC resolved, "
                        f"{n_mixed} remain ambiguous (AF {1-args.af_threshold:.2f}~{args.af_threshold:.2f})",
                        file=sys.stderr,
                    )

    # ── Splice exons ──────────────────────────────────────────────────────────
    spliced_nt = splice_exons(seq, seq_start0, exons, strand)

    # ── Translate ─────────────────────────────────────────────────────────────
    protein = translate(spliced_nt).rstrip('*')

    # ── Write output ──────────────────────────────────────────────────────────
    with open(args.out_nt, 'w') as fh:
        fh.write(f">{args.sample}|{args.gene}\n{spliced_nt}\n")

    with open(args.out_aa, 'w') as fh:
        fh.write(f">{args.sample}|{args.gene}\n{protein}\n")

    print(
        f"  {args.gene} [{args.sample}]: "
        f"{len(spliced_nt)} bp → {len(protein)} aa",
        file=sys.stderr,
    )

    # ── AF 표 저장 (--out-af 지정 시) ──────────────────────────────────────────
    if args.out_af is not None:
        write_af_table(af_records, args.out_af)


if __name__ == '__main__':
    main()
