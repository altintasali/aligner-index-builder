import os

import yaml
from snakemake.exceptions import WorkflowError
from snakemake.logging import logger
from snakemake.utils import validate

# -----------------------------------------------------------------------------
# Load & validate config
# -----------------------------------------------------------------------------
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
# Output paths (mirrors snakePipes createIndices' layout)
# -----------------------------------------------------------------------------
GENOME_FASTA = os.path.join(OUTDIR, "genome_fasta", "genome.fa")
GENOME_FAI = GENOME_FASTA + ".fai"
GENOME_2BIT = os.path.join(OUTDIR, "genome_fasta", "genome.2bit")
GENES_GTF = os.path.join(OUTDIR, "annotation", "genes.gtf") if GTF else None
TRANSCRIPTS_FA = os.path.join(OUTDIR, "annotation", "transcripts.fa")
CHROM_CHECK_OK = os.path.join(OUTDIR, "pipeline_info", "chrom_consistency.ok")

TOOL_TARGETS = {
    "star": os.path.join(OUTDIR, "STARIndex"),
    "hisat2": os.path.join(OUTDIR, "HISAT2Index", "genome.6.ht2"),
    "bowtie2": os.path.join(OUTDIR, "BowtieIndex", "genome.rev.2.bt2"),
    "bwa": os.path.join(OUTDIR, "BWAIndex", "genome.fa.sa"),
    "bwa-mem2": os.path.join(OUTDIR, "BWA-MEM2Index", "genome.fa.bwt.2bit.64"),
    "bwameth": os.path.join(OUTDIR, "BWAmethIndex", "genome.fa.bwameth.c2t.sa"),
    "bwameth2": os.path.join(OUTDIR, "BWAmeth2Index", "genome.fa.bwameth.c2t.bwt.2bit.64"),
    "salmon": os.path.join(OUTDIR, "SalmonIndex", "seq.bin"),
    # star and bismark rules declare directory() outputs, so the targets are
    # the directories themselves; the marker files asserted in CI/tests live
    # inside them (STARIndex/SAindex, .../CT_conversion/genome.1.bt2).
    "bismark": os.path.join(OUTDIR, "BismarkIndex", "Bismark_Genome"),
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
# -----------------------------------------------------------------------------
GENERATED_ENV_DIR = os.path.abspath("workflow/envs/generated")
os.makedirs(GENERATED_ENV_DIR, exist_ok=True)


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


STAR_ENV = _write_env("star", [f"star={V['star']}"])
SAMTOOLS_ENV = _write_env("samtools", [f"samtools={V['samtools']}"])
HISAT2_ENV = _write_env("hisat2", [f"hisat2={V['hisat2']}"])
BOWTIE2_ENV = _write_env("bowtie2", [f"bowtie2={V['bowtie2']}"])
BWA_ENV = _write_env("bwa", [f"bwa={V['bwa']}"])
BWA_MEM2_ENV = _write_env("bwa_mem2", [f"bwa-mem2={V['bwa_mem2']}"])
BWAMETH_ENV = _write_env(
    "bwa_meth",
    [f"bwa={V['bwa']}"]
    + ([f"bwa-meth={V['bwa_meth']}"] if V["bwa_meth"] else ["bwa-meth"]),
)
BISMARK_ENV = _write_env(
    "bismark", [f"bismark={V['bismark']}", f"bowtie2={V['bowtie2']}"]
)
GFFREAD_ENV = _write_env("gffread", [f"gffread={V['gffread']}"])
SALMON_ENV = _write_env("salmon", [f"salmon={V['salmon']}"])
UCSC_ENV = _write_env("ucsc", [f"ucsc-fatotwobit={V['ucsc_fatotwobit']}"])
PYTHON_ENV = _write_env("python", ["python>=3.9"])

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
