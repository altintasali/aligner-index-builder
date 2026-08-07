#!/usr/bin/env python3
"""Generate a tiny synthetic reference genome + GTF for CI / local smoke tests.

Deterministic (fixed seed). Writes:
  <outdir>/genome.fa   -- 3 contigs, ~120 kb of random DNA
  <outdir>/genome.gtf  -- a handful of genes (1 transcript, 2-4 exons each)

The dataset is deliberately small so the whole index set builds in a couple of
minutes on a laptop or in CI. STAR needs --genomeSAindexNbases /
--genomeChrBinNbits tuned down for references this small (its own formulas in
the manual, section 2.2.5); this script prints the recommended flags:

  genomeSAindexNbases = min(14, log2(genomeLength)/2 - 1)
  genomeChrBinNbits   = min(18, log2(genomeLength))

Usage:
  python .tests/generate_test_data.py [-o OUTDIR]
"""
import argparse
import math
import os
import random

CONTIGS = {"chr1": 50000, "chr2": 40000, "chr3": 30000}
SEED = 42
N_GENES = 60


def rand_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def write_fasta(path, contigs):
    with open(path, "w") as fh:
        for name, seq in contigs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i : i + 80] + "\n")


def write_gtf(path, genes):
    with open(path, "w") as fh:
        for g in genes:
            attrs = (
                f'gene_id "{g["gene_id"]}"; transcript_id '
                f'"{g["transcript_id"]}"; gene_name "{g["gene_name"]}";'
            )
            fh.write(
                f'{g["chrom"]}\ttest\tgene\t{g["start"]}\t{g["end"]}\t.\t'
                f'{g["strand"]}\t.\t{attrs}\n'
            )
            fh.write(
                f'{g["chrom"]}\ttest\ttranscript\t{g["start"]}\t{g["end"]}\t.\t'
                f'{g["strand"]}\t.\t{attrs}\n'
            )
            for start, end in g["exons"]:
                fh.write(
                    f'{g["chrom"]}\ttest\texon\t{start}\t{end}\t.\t'
                    f'{g["strand"]}\t.\t{attrs}\n'
                )


def make_gene(i, chrom, chrom_len, rng):
    strand = rng.choice(["+", "-"])
    n_exons = rng.randint(2, 4)
    exon_lens = [rng.randint(50, 150) for _ in range(n_exons)]
    intron_lens = [rng.randint(20, 80) for _ in range(n_exons - 1)]
    total = sum(exon_lens) + sum(intron_lens)
    start = rng.randint(1, max(1, chrom_len - total - 1))
    pos = start
    exons = []
    for j, el in enumerate(exon_lens):
        exons.append((pos, pos + el - 1))
        pos += el
        if j < len(intron_lens):
            pos += intron_lens[j]
    end = exons[-1][1]
    return {
        "chrom": chrom,
        "start": start,
        "end": end,
        "strand": strand,
        "exons": exons,
        "gene_id": f"GENE{i:04d}",
        "transcript_id": f"TX{i:04d}",
        "gene_name": f"gene_{i:04d}",
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-o", "--outdir", default=".tests/resources")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    rng = random.Random(SEED)

    contigs = {name: rand_seq(rng, length) for name, length in CONTIGS.items()}
    write_fasta(os.path.join(args.outdir, "genome.fa"), contigs)

    genes = []
    for i in range(N_GENES):
        chrom = rng.choice(sorted(CONTIGS))
        genes.append(make_gene(i, chrom, CONTIGS[chrom], rng))
    write_gtf(os.path.join(args.outdir, "genome.gtf"), genes)

    genome_len = sum(CONTIGS.values())
    sa_index = min(14, math.log2(genome_len) / 2 - 1)
    chr_bin = min(18, math.log2(genome_len))
    print(f"wrote {args.outdir}/genome.fa ({genome_len} bp) and genome.gtf")
    print(f"recommended STAR small-genome flags (for this fixture):")
    print(f"  --star-extra '--genomeSAindexNbases {int(sa_index)} "
          f"--genomeChrBinNbits {int(chr_bin)}'")


if __name__ == "__main__":
    main()
