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


def _masking_status(path):
    """Analyse genome masking from a FASTA file.

    Returns (status, soft_count, soft_pct, hard_count, hard_pct) where:
      status   = "unmasked" | "soft-masked" | "hard-masked" | "both"
      soft_count / hard_count = base counts
      soft_pct / hard_pct    = percentage of total bases
    """
    counts = {"upper": 0, "lower": 0, "n": 0, "other": 0}
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            for ch in line.strip():
                if ch in "ACGT":
                    counts["upper"] += 1
                elif ch in "acgt":
                    counts["lower"] += 1
                elif ch in "Nn":
                    counts["n"] += 1
                else:
                    counts["other"] += 1
    total = counts["upper"] + counts["lower"] + counts["n"] + counts["other"]
    if total == 0:
        return "unmasked", 0, 0.0, 0, 0.0
    soft_count = counts["lower"]
    hard_count = counts["n"]
    soft_pct = 100.0 * soft_count / total
    hard_pct = 100.0 * hard_count / total
    if soft_count > 0 and hard_count > 0:
        status = "both"
    elif soft_count > 0:
        status = "soft-masked"
    elif hard_count > 0:
        status = "hard-masked"
    else:
        status = "unmasked"
    return status, soft_count, soft_pct, hard_count, hard_pct



def _parse_gtf_features(path):
    """Count features by type (column 3) and genes per chromosome (column 1).

    Returns (feature_counts, chrom_gene_counts) where:
      feature_counts   = {feature_type: count}
      chrom_gene_counts = {chrom: gene_count}
    """
    from collections import Counter

    feature_counts = Counter()
    chrom_genes = Counter()
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 3:
                continue
            feature_counts[cols[2]] += 1
            if cols[2] == "gene":
                chrom_genes[cols[0]] += 1
    return dict(feature_counts), dict(chrom_genes)


def _count_gtf_features(path, feature):
    """Count occurrences of a specific feature type in a GTF."""
    feature_counts, _ = _parse_gtf_features(path)
    return feature_counts.get(feature, 0)


def _parse_transcript_lengths(gtf_path, fasta_path):
    """Compute transcript lengths by summing exon spans.

    Returns a dict of {transcript_id: length_bp} for every transcript in the
    GTF.  Transcripts without exon features are skipped.
    """
    from collections import defaultdict

    exons = defaultdict(list)
    with open(gtf_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9 or cols[2] != "exon":
                continue
            # Parse transcript_id from the 9th column attributes.
            attrs = cols[8]
            tx_id = None
            for tok in attrs.split(";"):
                tok = tok.strip()
                if tok.startswith("transcript_id"):
                    # transcript_id "ENST00000..."
                    tx_id = tok.split('"')[1] if '"' in tok else tok.split()[-1]
                    break
            if tx_id:
                exons[tx_id].append((int(cols[3]), int(cols[4])))
    # Length = sum of (end - start + 1) for each exon.
    return {tx: sum(e - s + 1 for s, e in spans) for tx, spans in exons.items()}


_TX_LENGTH_BINS = [
    ("0-500 bp", 0, 500),
    ("500-1 kb", 500, 1000),
    ("1-2 kb", 1000, 2000),
    ("2-5 kb", 2000, 5000),
    ("5-10 kb", 5000, 10000),
    ("10+ kb", 10000, float("inf")),
]


def _bin_transcript_lengths(lengths):
    """Bin transcript lengths into fixed categories.

    Returns an ordered {label: count} dict.
    """
    counts = {label: 0 for label, _, _ in _TX_LENGTH_BINS}
    for length in lengths.values():
        for label, lo, hi in _TX_LENGTH_BINS:
            if lo <= length < hi:
                counts[label] += 1
                break
    return counts


def _write_bargraph(path, doc_id, section_name, description, data, pconfig=None):
    """Write a MultiQC bargraph custom-content JSON.

    data: {sample_name: {category: count}}
    """
    if pconfig is None:
        pconfig = {}
    doc = {
        "id": doc_id,
        "section_name": section_name,
        "description": description,
        "plot_type": "bargraph",
        "pconfig": {
            "id": f"{doc_id}_plot",
            "title": section_name,
            "ylab": "Count",
            "cpswitch_counts_label": "Count",
            "cpswitch_catted_label": "Categorical",
            **pconfig,
        },
        "data": data,
    }
    with open(path, "w") as fh:
        json.dump(doc, fh, indent=2)


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
    import sys
    sys.stderr = open(snakemake.log[0], "w")  # noqa: F821

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
    }

    # Masking status -- always computed on the (filtered) genome; when
    # chromosome filtering was applied we also report the original.
    filtering_applied = params.get("filtering_applied", False)
    status, sc, sp, hc, hp = _masking_status(inp.fasta)
    if filtering_applied and not _is_empty(inp.raw_fasta):
        raw_status, raw_sc, raw_sp, raw_hc, raw_hp = _masking_status(
            inp.raw_fasta
        )
        resources.update({
            "Masking status (original)": raw_status,
            "Masking status (filtered)": status,
            "Soft-masked bases (original)": f"{raw_sc:,} ({raw_sp:.1f}%)",
            "Soft-masked bases (filtered)": f"{sc:,} ({sp:.1f}%)",
            "Hard-masked bases (N, original)": f"{raw_hc:,} ({raw_hp:.1f}%)",
            "Hard-masked bases (N, filtered)": f"{hc:,} ({hp:.1f}%)",
        })
    else:
        resources.update({
            "Masking status": status,
            "Soft-masked bases": f"{sc:,} ({sp:.1f}%)",
            "Hard-masked bases (N)": f"{hc:,} ({hp:.1f}%)",
        })

    resources.update({
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
    })
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
        _write_table(
            out.annotation,
            "annotation_summary",
            "Annotation summary",
            "Feature counts from the gene annotation used to build the "
            "splice-aware / transcriptome-aware indices.",
            annotation,
            col_header="Feature",
        )

    # GTF annotation plots (only when a GTF was provided).
    if has_gtf:
        feature_counts, chrom_gene_counts = _parse_gtf_features(inp.gtf)
        genome_label = params.get("genome_name", "genome")

        # Feature type bar chart.
        _write_bargraph(
            out.gtf_features,
            "gtf_features",
            "GTF feature types",
            "Count of each feature type in the gene annotation GTF.",
            {genome_label: feature_counts},
        )

        # Genes per chromosome bar chart.
        _write_bargraph(
            out.gtf_chroms,
            "gtf_chroms",
            "Genes per chromosome",
            "Number of gene features per chromosome/contig in the GTF.",
            {genome_label: chrom_gene_counts},
            pconfig={"cpswitch_counts_label": "Genes"},
        )

        # Transcript length distribution.
        tx_lengths = _parse_transcript_lengths(inp.gtf, inp.fasta)
        tx_bins = _bin_transcript_lengths(tx_lengths) if tx_lengths else {}
        _write_bargraph(
            out.gtf_tx_length,
            "gtf_tx_length",
            "Transcript length distribution",
            "Distribution of transcript lengths (sum of exon spans) from "
            "the GTF, binned into fixed size categories.",
            {genome_label: tx_bins} if tx_bins else {},
            pconfig={"cpswitch_counts_label": "Transcripts"},
        )


if __name__ == "__main__":
    main()
