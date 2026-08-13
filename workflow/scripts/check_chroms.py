"""Fail fast if the FASTA and GTF don't share any chromosome/contig name.

Runs as the check_chrom_consistency rule (see rules/index.smk): STAR, HISAT2
and salmon build their indices against the FASTA while using the GTF for
splice junctions/transcripts, so a mismatch in contig naming would silently
produce a broken index set. Exits non-zero (and the pipeline fails) unless at
least one contig name is shared. Also writes the per-side/shared contig counts
to {output.stats} (a JSON the report's resources_used rule reads).
"""
import json
import sys


def _fasta_contigs(path):
    contigs = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                contigs.add(line[1:].split()[0])
    return contigs


def _gtf_contigs(path):
    contigs = set()
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            contigs.add(line.split("\t", 1)[0])
    return contigs


def main():
    # snakemake's `script:` directive does not capture the process's output
    # into the rule log on its own -- redirect stderr to it so failure
    # messages are recorded where the user is told to look.
    sys.stderr = open(snakemake.log[0], "w")  # noqa: F821
    fasta_path = snakemake.input.fasta  # noqa: F821
    gtf_path = snakemake.input.gtf  # noqa: F821

    fasta = _fasta_contigs(fasta_path)
    gtf = _gtf_contigs(gtf_path)
    shared = fasta & gtf

    with open(snakemake.output.stats, "w") as fh:  # noqa: F821
        json.dump(
            {
                "fasta_contigs": len(fasta),
                "gtf_contigs": len(gtf),
                "shared_contigs": len(shared),
            },
            fh,
            indent=2,
        )

    if not shared:
        print(
            "check_chrom_consistency: FAILED -- no chromosome/contig name is "
            f"shared between the FASTA ({len(fasta)} contigs, e.g. "
            f"{sorted(fasta)[:5]}...) and the GTF ({len(gtf)} contigs, e.g. "
            f"{sorted(gtf)[:5]}...). The reference and annotation disagree on "
            "contig naming; refusing to build indices against them.",
            file=sys.stderr,
            flush=True,
        )
        sys.exit(1)

    print(
        f"check_chrom_consistency: OK -- {len(shared)} shared contig(s) "
        f"(e.g. {sorted(shared)[:5]}...).",
        file=sys.stderr,
        flush=True,
    )


if __name__ == "__main__":
    main()
