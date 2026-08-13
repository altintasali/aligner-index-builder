import os

import snakemake
import yaml
from snakemake.exceptions import WorkflowError
from snakemake.logging import logger
from snakemake.utils import validate

# Repository layout anchors. `repo_root` comes from the run config (written by
# the CLI), so nothing depends on the directory snakemake is invoked from --
# the CLI runs snakemake with `--directory <outdir>` so the `.snakemake`
# metadata/lock dir lands in the output dir, never the repo. Direct
# hand-written configs fall back to the cwd (run them from the repo root).
REPO_ROOT = config.get("repo_root") or os.getcwd()
WORKFLOW_DIR = os.path.join(REPO_ROOT, "workflow")

# -----------------------------------------------------------------------------
# Load & validate config
# -----------------------------------------------------------------------------
# Tool names are case-insensitive: normalize to lowercase before validation so
# direct configs accept `tools: [STAR]` or `tools: [All]` just like the CLI
# does. Non-string entries are left alone for the schema to reject.
_raw_tools = config.get("tools")
if isinstance(_raw_tools, str):
    config["tools"] = _raw_tools.lower()
elif isinstance(_raw_tools, list):
    config["tools"] = [t.lower() if isinstance(t, str) else t for t in _raw_tools]

validate(config, schema="../schemas/config.schema.yaml")

V = config["versions"]

OUTDIR = config["outdir"]
GENOME_NAME = config.get("genome_name", "genome")


def _is_url(path):
    return "://" in str(path)


FASTA = config["fasta"]
GTF = config.get("gtf") or None

# Reference inputs may be local paths (plain or gzipped) or URLs. Local paths
# are resolved against the directory snakemake is run from (the CLI resolves
# them to absolute paths before writing the config); fail fast at parse time
# on a missing/typo'd path instead of letting a job die hours into a queue.
for _key, _path in (("fasta", FASTA), ("gtf", GTF)):
    if _path and not _is_url(_path) and not os.path.exists(str(_path)):
        raise WorkflowError(
            f"Missing reference file ({_key}): {_path!r} does not exist "
            f"(resolves to: {os.path.abspath(str(_path))}).\n"
            "Fix the path in the config, or pass a URL. The CLI resolves "
            "local paths to absolute paths automatically."
        )

# -----------------------------------------------------------------------------
# Tools: which indices to build. `all` is expanded by the CLI; a direct
# config may list tools explicitly. Tools that need the GTF (STAR, HISAT2,
# salmon) are dropped with a warning when no GTF is given -- the DNA aligners
# (bowtie2, bwa, bwa-mem2, bwameth/bwameth2, bismark) only need the FASTA.
# -----------------------------------------------------------------------------
ALL_TOOLS = [
    "star", "hisat2", "bowtie2", "bwa", "bwa-mem2",
    "bwameth", "bwameth2", "salmon", "bismark",
]
GTF_TOOLS = {"star", "hisat2", "salmon"}

_tools = config["tools"]
if _tools == "all":
    _tools = list(ALL_TOOLS)
_tools = [t for t in _tools if t in ALL_TOOLS]

if not GTF:
    _dropped = [t for t in _tools if t in GTF_TOOLS]
    if _dropped:
        logger.warning(
            "No GTF given; skipping GTF-dependent tools: "
            + ", ".join(_dropped)
        )
        _tools = [t for t in _tools if t not in GTF_TOOLS]

if not _tools:
    raise WorkflowError(
        "No tools left to build (all requested tools need a GTF, but none "
        "was given)."
    )

TOOLS = _tools

# -----------------------------------------------------------------------------
# Output paths (all per-tool index dirs live under <outdir>/index/<tool>/)
# -----------------------------------------------------------------------------
INDEX_DIR = os.path.join(OUTDIR, "index")
GENOME_FASTA = os.path.join(OUTDIR, "genome_fasta", "genome.fa")
GENOME_FAI = GENOME_FASTA + ".fai"
GENOME_2BIT = os.path.join(OUTDIR, "genome_fasta", "genome.2bit")
GENES_GTF = os.path.join(OUTDIR, "annotation", "genes.gtf") if GTF else None
TRANSCRIPTS_FA = os.path.join(OUTDIR, "annotation", "transcripts.fa")
CHROM_CHECK_OK = os.path.join(OUTDIR, "pipeline_info", "chrom_consistency.ok")
CHROM_STATS = os.path.join(OUTDIR, "pipeline_info", "chrom_consistency.json")

TOOL_TARGETS = {
    "star": os.path.join(INDEX_DIR, "star"),
    "hisat2": os.path.join(INDEX_DIR, "hisat2", "genome.6.ht2"),
    "bowtie2": os.path.join(INDEX_DIR, "bowtie2", "genome.rev.2.bt2"),
    "bwa": os.path.join(INDEX_DIR, "bwa", "genome.fa.sa"),
    "bwa-mem2": os.path.join(INDEX_DIR, "bwa-mem2", "genome.fa.bwt.2bit.64"),
    "bwameth": os.path.join(INDEX_DIR, "bwameth", "genome.fa.bwameth.c2t.sa"),
    "bwameth2": os.path.join(INDEX_DIR, "bwameth2", "genome.fa.bwameth.c2t.bwt.2bit.64"),
    "salmon": os.path.join(INDEX_DIR, "salmon", "seq.bin"),
    # star and bismark rules declare directory() outputs, so the targets are
    # the directories themselves; the marker files asserted in CI/tests live
    # inside them (index/star/SAindex, index/bismark/Bisulfite_Genome/...).
    "bismark": os.path.join(INDEX_DIR, "bismark", "Bisulfite_Genome"),
}


def all_index_targets():
    """The complete list of outputs the run should produce."""
    targets = [GENOME_FASTA, GENOME_FAI, GENOME_2BIT]
    if GENES_GTF:
        targets.append(CHROM_CHECK_OK)
    for tool in TOOLS:
        targets.append(TOOL_TARGETS[tool])
    return targets


# -----------------------------------------------------------------------------
# Per-rule compute resources, from workflow/default-config/resources.yaml
# (optionally overridden by a `resources:` block in the run's config).
# -----------------------------------------------------------------------------
RESOURCES = config.get("resources", {})
_RESOURCE_DEFAULTS = {"threads": 1, "mem_mb": 4000, "runtime": 60}


def get_resources(rule_name):
    return {**_RESOURCE_DEFAULTS, **RESOURCES.get(rule_name, {})}

# -----------------------------------------------------------------------------
# Generated per-tool conda env files (used with `--use-conda` / `--sdm
# conda`). The pre-built env ships the same tools, so plain `snakemake` runs
# never touch these -- they exist so a direct snakemake run can resolve its
# own tools if you prefer not to use the pre-built env.
# Absolute path: rules included from workflow/rules/ resolve relative paths
# against that directory, not the run directory.
GENERATED_ENV_DIR = os.path.join(WORKFLOW_DIR, "envs", "generated")
os.makedirs(GENERATED_ENV_DIR, exist_ok=True)


def _pin(name, version):
    """A conda spec `name=version` (or bare `name` when unpinned). Versions
    that already carry an operator (e.g. gffread `>=0.12.6`) are appended
    verbatim; naive `name={version}` would render `gffread=>=0.12.6`, which
    conda rejects."""
    if not version:
        return name
    if version[0] in "=<>":
        return f"{name}{version}"
    return f"{name}={version}"


def _write_env(name, dependencies):
    path = f"{GENERATED_ENV_DIR}/{name}.yaml"
    with open(path, "w") as fh:
        yaml.safe_dump(
            {
                "channels": ["bioconda", "conda-forge"],
                "dependencies": list(dependencies),
            },
            fh,
            sort_keys=False,
        )
    return path


STAR_ENV = _write_env("star", [_pin("star", V["star"])])
SAMTOOLS_ENV = _write_env("samtools", [_pin("samtools", V["samtools"])])
HISAT2_ENV = _write_env("hisat2", [_pin("hisat2", V["hisat2"])])
BOWTIE2_ENV = _write_env("bowtie2", [_pin("bowtie2", V["bowtie2"])])
BWA_ENV = _write_env("bwa", [_pin("bwa", V["bwa"])])
BWA_MEM2_ENV = _write_env("bwa_mem2", [_pin("bwa-mem2", V["bwa_mem2"])])
BWAMETH_ENV = _write_env(
    "bwa_meth",
    [_pin("bwa", V["bwa"]), _pin("bwa-mem2", V["bwa_mem2"])]
    + [_pin("bwameth", V["bwameth"])],
)
BISMARK_ENV = _write_env(
    "bismark", [_pin("bismark", V["bismark"]), _pin("bowtie2", V["bowtie2"])]
)
GFFREAD_ENV = _write_env("gffread", [_pin("gffread", V["gffread"])])
SALMON_ENV = _write_env("salmon", [_pin("salmon", V["salmon"])])
UCSC_ENV = _write_env("ucsc", [_pin("ucsc-fatotwobit", V["ucsc_fatotwobit"])])
PYTHON_ENV = _write_env("python", ["python>=3.9"])
MULTIQC_ENV = _write_env("multiqc", [_pin("multiqc", V["multiqc"])])

# -----------------------------------------------------------------------------
# Versions recorded in the MultiQC report / Software Versions section
# -----------------------------------------------------------------------------
# Pipeline version: read from the VERSION file at the repo root (same value
# the CLI echoes with --version).
try:
    with open(os.path.join(REPO_ROOT, "VERSION")) as _fh:
        PIPELINE_VERSION = _fh.read().strip()
except OSError:
    PIPELINE_VERSION = "unknown"

SNAKEMAKE_VERSION = snakemake.__version__

# Index directories, keyed by tool name -- used by the resources_used rule to
# report per-index disk usage. Mirrors the CLI's TOOL_DIRS mapping.
TOOL_INDEX_DIRS = {
    "star": os.path.join(INDEX_DIR, "star"),
    "hisat2": os.path.join(INDEX_DIR, "hisat2"),
    "bowtie2": os.path.join(INDEX_DIR, "bowtie2"),
    "bwa": os.path.join(INDEX_DIR, "bwa"),
    "bwa-mem2": os.path.join(INDEX_DIR, "bwa-mem2"),
    "bwameth": os.path.join(INDEX_DIR, "bwameth"),
    "bwameth2": os.path.join(INDEX_DIR, "bwameth2"),
    "salmon": os.path.join(INDEX_DIR, "salmon"),
    "bismark": os.path.join(INDEX_DIR, "bismark"),
}

# -----------------------------------------------------------------------------
# MultiQC report outputs (qc.smk)
# -----------------------------------------------------------------------------
# Absolute path so the multiqc shell command works no matter what directory
# snakemake is invoked from.
MULTIQC_CONFIG = os.path.join(WORKFLOW_DIR, "default-config", "multiqc_config.yaml")

QC_VERSIONS = os.path.join(OUTDIR, "versions", f"{GENOME_NAME}_mqc_versions.yml")
QC_RESOURCES_USED = os.path.join(OUTDIR, "pipeline_info", "resources_used_mqc.json")
QC_CONFIG_USED = os.path.join(OUTDIR, "pipeline_info", "config_used_mqc.json")
QC_INDEX_SIZES = os.path.join(OUTDIR, "pipeline_info", "index_sizes_mqc.json")
QC_ANNOTATION_SUMMARY = os.path.join(
    OUTDIR, "pipeline_info", "annotation_summary_mqc.json"
)
QC_BENCHMARK_SUMMARY = os.path.join(
    OUTDIR, "pipeline_info", "benchmark_summary_mqc.json"
)
QC_MULTIQC_HTML = os.path.join(OUTDIR, "qc", "multiqc_report.html")
QC_MULTIQC_DATA = os.path.join(OUTDIR, "qc", "multiqc_report_data")


def all_benchmark_files():
    """Every benchmark file this run will produce, for the benchmark_summary
    rule (qc.smk) to aggregate into the MultiQC resource-usage section. Built
    deterministically from TOOLS/GTF (rules only run when their tool is
    requested); the report rules' own benchmarks (software_versions,
    resources_used, config_used, benchmark_summary, multiqc) are excluded --
    they are negligible and would otherwise create a dependency cycle."""
    _qc_rules = {
        "software_versions", "resources_used", "config_used",
        "benchmark_summary", "multiqc",
    }
    _dir = os.path.join(OUTDIR, "pipeline_info", "benchmarks")
    _to_rule = {
        "star": "star",
        "bowtie2": "bowtie2",
        "bwa": "bwa",
        "bwa-mem2": "bwa_mem2",
        "bwameth": "bwameth",
        "bwameth2": "bwameth2",
        "bismark": "bismark",
    }
    files = [
        "prepare_genome", "fasta_fai", "fasta_2bit",
    ]
    if GTF:
        files += ["prepare_gtf", "chrom_consistency"]
    for tool in TOOLS:
        if tool == "hisat2":
            files += ["hisat2_splicesites", "hisat2_exons", "hisat2"]
        elif tool == "salmon":
            files += ["transcriptome_fasta", "salmon"]
        else:
            files.append(_to_rule[tool])
    return sorted(
        os.path.join(_dir, stem + ".txt")
        for stem in files
        if stem not in _qc_rules
    )


def allocated_resources_by_rule():
    """{rule: {"threads", "mem_mb"}} for every rule with a benchmark file,
    read from resources.yaml -- the per-job allocation against which the
    benchmark_summary script computes CPU/RAM efficiency. Iterated in
    Snakemake's own rule order so the resource-usage table lists rules in
    workflow order, not alphabetically."""
    _benchmark_rules = {
        os.path.basename(path)[:-4] for path in all_benchmark_files()
    }
    _out = {}
    for _r in workflow.rules:
        if _r.name in _benchmark_rules:
            _out[_r.name] = get_resources(_r.name)
    return _out

# -----------------------------------------------------------------------------
# Lightweight record of how this run's config resolved, for the user to check
# -----------------------------------------------------------------------------
os.makedirs(os.path.join(OUTDIR, "pipeline_info", "logs"), exist_ok=True)
with open(os.path.join(OUTDIR, "pipeline_info", "logs", "config_resolution.log"), "w") as fh:
    fh.write(f"genome_name = {GENOME_NAME}\n")
    fh.write(f"outdir = {OUTDIR}\n")
    fh.write(f"fasta = {FASTA}\n")
    fh.write(f"gtf = {GTF if GTF else '(none)'}\n")
    fh.write(f"tools = {', '.join(TOOLS)}\n")
    fh.write(f"sjdb_overhang = {config.get('sjdb_overhang', 100)}\n")
    fh.write(f"salmon_kmer = {config.get('salmon_kmer', 31)}\n")
