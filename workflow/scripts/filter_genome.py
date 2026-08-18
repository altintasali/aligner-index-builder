"""Filter contigs from a FASTA based on chromosome name rules.

Runs as the filter_genome rule (see rules/index.smk): optionally removes
non-standard chromosomes and/or applies user-supplied keep/exclude regex
patterns before any index is built.  When no filtering is requested the
rule is a trivial copy (the DAG always includes it so downstream rules
read the same GENOME_FASTA path regardless).
"""
import re
import sys

# Built-in pattern for "standard" chromosomes -- covers UCSC (chr prefix,
# lowercase) and Ensembl (no prefix, uppercase) naming for human, mouse,
# pig, chicken and most model organisms.  case-insensitive to handle
# edge-case assemblies that mix conventions (e.g. "Chr1", "CHRMT").
_STANDARD_RE = re.compile(r"^(chr)?(\d+|[XYWZ]|MT?)$", re.IGNORECASE)


def _compile(patterns):
    """Compile a list of regex strings into compiled pattern objects."""
    return [re.compile(p) for p in patterns]


def _keep(contig, keep_res, remove_standard, exclude_res):
    """Return True if *contig* should be kept."""
    # 1. --keep-chroms: positive selection (applied first)
    if keep_res and not any(p.search(contig) for p in keep_res):
        return False
    # 2. --remove-nonstandard: built-in filter
    if remove_standard and not _STANDARD_RE.match(contig):
        return False
    # 3. --exclude-chroms: custom exclusion
    if exclude_res and any(p.search(contig) for p in exclude_res):
        return False
    return True


def _fasta_iter(path):
    """Yield (header_line, [sequence_lines]) tuples from a FASTA file."""
    header = None
    seq_lines = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, seq_lines
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            yield header, seq_lines


def main():
    sys.stderr = open(snakemake.log[0], "w")  # noqa: F821
    inp = str(snakemake.input)  # noqa: F821
    out = str(snakemake.output)  # noqa: F821
    params = snakemake.params  # noqa: F821

    remove_standard = params.get("remove_nonstandard", False)
    keep_patterns = params.get("keep_chroms", [])
    exclude_patterns = params.get("exclude_chroms", [])

    keep_res = _compile(keep_patterns) if keep_patterns else []
    exclude_res = _compile(exclude_patterns) if exclude_patterns else []

    if not keep_res and not remove_standard and not exclude_res:
        # No filtering requested -- trivial copy
        import shutil
        shutil.copy2(inp, out)
        print(
            "filter_genome: no filtering requested -- copied input to output",
            file=sys.stderr,
            flush=True,
        )
        return

    kept = 0
    removed = 0
    removed_names = []
    with open(out, "w") as fout:
        for header, seq_lines in _fasta_iter(inp):
            contig = header[1:].split()[0]
            if _keep(contig, keep_res, remove_standard, exclude_res):
                fout.write(header)
                fout.writelines(seq_lines)
                kept += 1
            else:
                removed += 1
                if removed <= 20:
                    removed_names.append(contig)

    summary = (
        f"filter_genome: kept {kept} contig(s), removed {removed} contig(s)"
    )
    if removed_names:
        preview = ", ".join(removed_names)
        if removed > 20:
            preview += ", ..."
        summary += f" (e.g. {preview})"
    print(summary, file=sys.stderr, flush=True)

    if kept == 0:
        print(
            "filter_genome: ERROR -- no contigs survived filtering. "
            "Check your --keep-chroms / --remove-nonstandard / --exclude-chroms "
            "arguments.",
            file=sys.stderr,
            flush=True,
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
