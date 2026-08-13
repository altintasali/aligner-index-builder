#!/usr/bin/env python3
"""Write the "Resources used", "Index sizes" and "Annotation summary" custom
content sections of the MultiQC report (qc.smk -> resources_used rule).

Everything derives from files the run already produced:
  - SHA-256 of the normalized reference (genome.fa) and annotation (genes.gtf)
    for reproducibility;
  - contig count / total length / N50 from the .fai and GC% from a single
    pass over genome.fa;
  - the FASTA-vs-GTF contig consistency result from chrom_consistency.json
    (written by check_chroms.py);
  - per-index disk usage (os.walk over each index dir; symlinks count as the
    link itself so the shared genome.fa symlink is not double-counted);
  - gene/transcript counts from the GTF and the gffread transcriptome.

Each output is a self-describing MultiQC custom content JSON (id +
section_name + plot_type: table). Section ids are referenced by
custom_content.order in multiqc_config.yaml.
"""
import hashlib
import json
import os

try:
    # Only defined when run by Snakemake's `script:` directive.
    snakemake
except NameError:  # pragma: no cover
    snakemake = None


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _fai_stats(path):
    """(count, total_bp, n50_bp, max_bp) from a samtools faidx .fai."""
    lengths = []
    with open(path) as fh:
        for line in fh:
            fields = line.split("\t")
            if fields:
                try:
                    lengths.append(int(fields[1]))
                except ValueError:
                    pass
    if not lengths:
        return 0, 0, 0, 0
    lengths.sort(reverse=True)
    total = sum(lengths)
    half = total / 2.0
    cum = 0
    n50 = lengths[-1]
    for length in lengths:
        cum += length
        if cum >= half:
            n50 = length
            break
    return len(lengths), total, n50, lengths[0]


def _gc_percent(path):
    counts = {"a": 0, "c": 0, "g": 0, "t": 0}
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            for ch in line.strip().lower():
                if ch in counts:
                    counts[ch] += 1
    denom = sum(counts.values())
    if not denom:
        return 0.0
    return 100.0 * (counts["g"] + counts["c"]) / denom


def _count_gtf_features(path, feature):
    n = 0
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            if line.split("\t")[2] == feature:
                n += 1
    return n


def _count_headers(path):
    n = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def _dir_size(path):
    total = 0
    for root, _dirs, files in os.walk(path):
        for name in files:
            try:
                total += os.lstat(os.path.join(root, name)).st_size
            except OSError:
                pass
    return total


def _human(n):
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if n < 1024:
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} PiB"


def _write_table(path, doc_id, section_name, description, rows, col_header="Resource"):
    """rows: ordered {row_name: {column_key: value}} with a shared headers map."""
    doc = {
        "id": doc_id,
        "section_name": section_name,
        "description": description,
        "plot_type": "table",
        "pconfig": {
            "id": f"{doc_id}_table",
            "title": section_name,
            "col1_header": col_header,
            "sort_rows": False,
        },
        "headers": {
            "value": {
                "title": "Value",
            },
        },
        "data": {k: {"value": v} for k, v in rows.items()},
    }
    with open(path, "w") as fh:
        json.dump(doc, fh, indent=2)


def _is_empty(value):
    return isinstance(value, list) and not value


def main():
    params = snakemake.params  # noqa: F821
    inp = snakemake.input  # noqa: F821
    out = snakemake.output  # noqa: F821

    fasta_path = params.fasta_path or inp.fasta
    has_gtf = not _is_empty(inp.gtf)
    gtf_path = params.gtf_path or ("(none)" if not has_gtf else inp.gtf)

    fasta_sha = _sha256(inp.fasta)
    n_contigs, total_bp, n50, max_bp = _fai_stats(inp.fai)
    gc = round(_gc_percent(inp.fasta), 3)

    index_sizes = {}
    for tool in params.tools:
        index_sizes[tool] = _dir_size(params.index_dirs[tool])
    total_index = sum(index_sizes.values())

    chrom = None
    if not _is_empty(inp.chrom):
        with open(inp.chrom) as fh:
            chrom = json.load(fh)

    resources = {
        "Reference genome (FASTA)": fasta_path,
        "Reference FASTA SHA-256": fasta_sha,
        "Genome contigs": n_contigs,
        "Genome total length (bp)": total_bp,
        "Genome N50 (bp)": n50,
        "Genome GC (%)": gc,
        "Gene annotation (GTF)": gtf_path if has_gtf else "(none)",
        "GTF SHA-256": _sha256(inp.gtf) if has_gtf else "(none)",
        "FASTA vs GTF contigs": (
            f"OK: {chrom['shared_contigs']} shared "
            f"(FASTA {chrom['fasta_contigs']} / GTF {chrom['gtf_contigs']})"
            if chrom else "not checked (no GTF)"
        ),
        "Index sets built": ", ".join(params.tools) or "(none)",
        "Total index size": _human(total_index),
        "Pipeline version": params.pipeline_version,
        "Snakemake version": params.snakemake_version,
    }
    _write_table(
        out.resources,
        "resources_used",
        "Resources used",
        "The reference genome and annotation this index set was built from, "
        "with checksums and sequence statistics, the index sets produced, and "
        "the pipeline/tool versions.",
        resources,
    )

    index_rows = {
        tool: _human(size) for tool, size in index_sizes.items()
    }
    index_rows["total"] = _human(total_index)
    _write_table(
        out.index_sizes,
        "index_sizes",
        "Index sizes",
        "Disk footprint of each built index directory (symlinked reference "
        "files count as the link itself).",
        index_rows,
        col_header="Index",
    )

    if not _is_empty(out.annotation):
        gtf = inp.gtf
        annotation = {
            "Genes (GTF)": _count_gtf_features(gtf, "gene"),
            "Transcripts (GTF)": _count_gtf_features(gtf, "transcript"),
        }
        if not _is_empty(inp.transcripts):
            annotation["Transcripts (gffread)"] = _count_headers(inp.transcripts)
        _write_table(
            out.annotation,
            "annotation_summary",
            "Annotation summary",
            "Feature counts from the gene annotation used to build the "
            "splice-aware / transcriptome-aware indices.",
            annotation,
            col_header="Feature",
        )


if __name__ == "__main__":
    main()
