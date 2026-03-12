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


if __name__ == '__main__':
    main()
