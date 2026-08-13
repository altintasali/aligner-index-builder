# =============================================================================
# MultiQC report rules. The pipeline has no per-sample data (no alignments or
# reads), so the report is entirely custom content plus the built-in Software
# Versions section:
#
#   software_versions  -> results/versions/<genome>_mqc_versions.yml
#                         (MultiQC's Software Versions section, auto-discovered)
#   resources_used     -> results/pipeline_info/{resources_used,index_sizes,
#                         annotation_summary}_mqc.json  ("Resources used",
#                         "Index sizes", "Annotation summary" tables)
#   config_used        -> results/pipeline_info/config_used_mqc.json
#                         ("Configuration used" table)
#   benchmark_summary  -> results/pipeline_info/benchmark_summary_mqc.json
#                         ("Resource usage" table, aggregated from the per-rule
#                         benchmark files)
#   multiqc            -> results/qc/multiqc_report.html (+ _data dir)
#
# The custom-content section order is set in
# workflow/default-config/multiqc_config.yaml.
# =============================================================================

# -----------------------------------------------------------------------------
# Software Versions section (MultiQC auto-discovers *_mqc_versions.yml)
# -----------------------------------------------------------------------------
rule software_versions:
    # Writes the pinned tool versions in the format MultiQC's Software
    # Versions section expects (a `*_mqc_versions.yml` file -- that exact
    # filename pattern, per MultiQC's search pattern). The versions come from
    # workflow/default-config/versions.yaml, the single source of truth.
    output:
        QC_VERSIONS,
    threads: get_resources("software_versions")["threads"]
    resources:
        mem_mb=get_resources("software_versions")["mem_mb"],
        runtime=get_resources("software_versions")["runtime"],
    benchmark:
        os.path.join(OUTDIR, "pipeline_info", "benchmarks", "software_versions.txt"),
    log:
        os.path.join(OUTDIR, "pipeline_info", "logs", "software_versions.log"),
    conda:
        PYTHON_ENV
    run:
        import yaml

        versions = {
            "STAR": V["star"],
            "SAMtools": V["samtools"],
            "HISAT2": V["hisat2"],
            "bowtie2": V["bowtie2"],
            "BWA": V["bwa"],
            "bwa-mem2": V["bwa_mem2"],
            "bwameth": V["bwameth"] if V["bwameth"] else "unpinned (latest)",
            "Bismark": V["bismark"],
            "salmon": V["salmon"],
            "gffread": V["gffread"],
            "faToTwoBit": V["ucsc_fatotwobit"],
            "MultiQC": V["multiqc"],
        }
        with open(QC_VERSIONS, "w") as fh:
            yaml.safe_dump(
                {"software_versions": versions}, fh, default_flow_style=False
            )


# -----------------------------------------------------------------------------
# "Resources used" / "Index sizes" / "Annotation summary" sections
# -----------------------------------------------------------------------------
rule resources_used:
    # Reference checksums + sequence statistics, the FASTA-vs-GTF consistency
    # result, per-index disk usage, and annotation feature counts -- everything
    # for the report's reference sections, written by resources_used.py. The
    # index marker files are inputs purely to order the rule after the builds;
    # params.index_dirs holds the actual directories to size.
    input:
        fasta=GENOME_FASTA,
        fai=GENOME_FAI,
        index_markers=[TOOL_TARGETS[t] for t in TOOLS],
        gtf=(GENES_GTF if GENES_GTF else []),
        chrom=(CHROM_STATS if GENES_GTF else []),
        transcripts=(TRANSCRIPTS_FA if "salmon" in TOOLS else []),
    output:
        resources=QC_RESOURCES_USED,
        index_sizes=QC_INDEX_SIZES,
        annotation=(QC_ANNOTATION_SUMMARY if GENES_GTF else []),
    params:
        fasta_path=FASTA,
        gtf_path=GTF or "",
        tools=TOOLS,
        index_dirs={t: TOOL_INDEX_DIRS[t] for t in TOOLS},
        pipeline_version=PIPELINE_VERSION,
        snakemake_version=SNAKEMAKE_VERSION,
    threads: get_resources("resources_used")["threads"]
    resources:
        mem_mb=get_resources("resources_used")["mem_mb"],
        runtime=get_resources("resources_used")["runtime"],
    benchmark:
        os.path.join(OUTDIR, "pipeline_info", "benchmarks", "resources_used.txt"),
    log:
        os.path.join(OUTDIR, "pipeline_info", "logs", "resources_used.log"),
    conda:
        PYTHON_ENV
    script:
        "../scripts/resources_used.py"


# -----------------------------------------------------------------------------
# "Configuration used" section
# -----------------------------------------------------------------------------
rule config_used:
    # The resolved run config, written by the CLI / read from --configfile,
    # so the report shows exactly what was requested.
    output:
        QC_CONFIG_USED,
    threads: get_resources("config_used")["threads"]
    resources:
        mem_mb=get_resources("config_used")["mem_mb"],
        runtime=get_resources("config_used")["runtime"],
    benchmark:
        os.path.join(OUTDIR, "pipeline_info", "benchmarks", "config_used.txt"),
    log:
        os.path.join(OUTDIR, "pipeline_info", "logs", "config_used.log"),
    conda:
        PYTHON_ENV
    run:
        import json

        rows = {
            "genome_name": GENOME_NAME,
            "outdir": OUTDIR,
            "fasta": FASTA,
            "gtf": GTF if GTF else "(none)",
            "tools": ", ".join(TOOLS),
            "sjdb_overhang": config.get("sjdb_overhang", 100),
            "star_extra": config.get("star_extra") or "(none)",
            "salmon_kmer": config.get("salmon_kmer", 31),
            "pipeline_version": PIPELINE_VERSION,
            "snakemake_version": SNAKEMAKE_VERSION,
        }
        doc = {
            "id": "config_used",
            "section_name": "Configuration used",
            "description": "The resolved run configuration (config file "
                           "generated by the aligner-index-builder CLI or "
                           "passed with --configfile).",
            "plot_type": "table",
            "pconfig": {
                "id": "config_used_table",
                "title": "Configuration used",
                "col1_header": "Setting",
                "sort_rows": False,
            },
            "headers": {"value": {"title": "Value"}},
            "data": {k: {"value": v} for k, v in rows.items()},
        }
        with open(QC_CONFIG_USED, "w") as fh:
            json.dump(doc, fh, indent=2)


# -----------------------------------------------------------------------------
# "Resource usage" section (aggregated per-rule benchmarks)
# -----------------------------------------------------------------------------
rule benchmark_summary:
    # Aggregates every benchmark file this run produces into a
    # self-describing MultiQC custom content table ("Resource usage",
    # id resource_usage), rendering the allocated-vs-used CPU/RAM picture.
    input:
        all_benchmark_files(),
    output:
        QC_BENCHMARK_SUMMARY,
    params:
        allocated=allocated_resources_by_rule(),
    threads: get_resources("benchmark_summary")["threads"]
    resources:
        mem_mb=get_resources("benchmark_summary")["mem_mb"],
        runtime=get_resources("benchmark_summary")["runtime"],
    benchmark:
        os.path.join(OUTDIR, "pipeline_info", "benchmarks", "benchmark_summary.txt"),
    log:
        os.path.join(OUTDIR, "pipeline_info", "logs", "benchmark_summary.log"),
    conda:
        PYTHON_ENV
    script:
        "../scripts/benchmark_summary.py"


# -----------------------------------------------------------------------------
# The report itself
# -----------------------------------------------------------------------------
rule multiqc:
    # Runs MultiQC over results/versions (the *_mqc_versions.yml) and
    # results/pipeline_info (the *_mqc.json custom content), applying the
    # custom config workflow/default-config/multiqc_config.yaml (title,
    # section order).
    input:
        versions=QC_VERSIONS,
        resources=QC_RESOURCES_USED,
        config=QC_CONFIG_USED,
        index_sizes=QC_INDEX_SIZES,
        annotation=(QC_ANNOTATION_SUMMARY if GENES_GTF else []),
        benchmark=QC_BENCHMARK_SUMMARY,
    output:
        html=QC_MULTIQC_HTML,
        data=directory(QC_MULTIQC_DATA),
    params:
        indirs=lambda wc, input: sorted({os.path.dirname(f) for f in input}),
    threads: get_resources("multiqc")["threads"]
    resources:
        mem_mb=get_resources("multiqc")["mem_mb"],
        runtime=get_resources("multiqc")["runtime"],
    benchmark:
        os.path.join(OUTDIR, "pipeline_info", "benchmarks", "multiqc.txt"),
    log:
        os.path.join(OUTDIR, "pipeline_info", "logs", "multiqc.log"),
    conda:
        MULTIQC_ENV
    shell:
        "multiqc {params.indirs} --force "
        "-c {MULTIQC_CONFIG} "
        "-o {OUTDIR}/qc -n multiqc_report "
        "> {log} 2>&1"
