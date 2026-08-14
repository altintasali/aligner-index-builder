# aligner-index-builder

A [Snakemake](https://snakemake.readthedocs.io/) pipeline that builds a
complete set of aligner/indexer indices from a reference genome **FASTA**
(+ optional **GTF**), inspired by the
[snakePipes](https://github.com/maxplanck-ie/snakepipes) `createIndices`
layout: STAR, HISAT2, bowtie2, bwa, bwa-mem2, bwameth, bwameth2,
**Bismark** and salmon -- plus the `.fai` / `.2bit`
side products and a `gffread` transcriptome for decoy-aware salmon indexes.
Each tool's index set lands in its own lowercase dir under `<outdir>/index/`
(`index/star/`, `index/bwa-mem2/`, `index/bwameth/`, ... -- see
[What gets built](#what-gets-built)).
Every run also produces a [MultiQC](https://multiqc.info/) report
(`qc/multiqc_report.html`) summarising the reference and annotation, the
per-index disk footprints, the resolved run configuration, per-rule resource
usage, and the pinned tool versions (see [MultiQC report](#multiqc-report)).

It is driven by a small CLI (`workflow/scripts/aligner-index-builder`) that
writes a run config and hands off to snakemake; tool versions are pinned
exactly. Everything runs out of a pre-built conda environment, distributed as
a GitHub Release asset -- no conda, no container runtime, no per-tool solves.

Companion of [TEtranscripts-pipe](https://github.com/altintasali/TEtranscripts-pipe)
(the RNA-seq quantification pipeline that consumes these index sets).

---

## Table of contents

- [Requirements](#requirements)
- [Install](#install)
- [Usage](#usage)
- [HPC / SLURM](#hpc--slurm)
- [What gets built](#what-gets-built)
- [Output layout](#output-layout)
- [MultiQC report](#multiqc-report)
- [Direct snakemake use](#direct-snakemake-use)
- [Config reference](#config-reference)
- [Architecture](#architecture)
- [Versioning & releases](#versioning--releases)
- [Caveats](#caveats)

## Requirements

- Linux **x86_64** for the pre-built environment (the release tarball is
  packed on GitHub Actions ubuntu runners). On other platforms create the
  environment yourself with conda (see below).
- The workflow itself needs only `bash`, `curl` and Python >= 3.9 for
  `install-env.sh` / the CLI.

## Install

### Option A -- pre-built environment (recommended)

```bash
git clone https://github.com/altintasali/aligner-index-builder.git
cd aligner-index-builder

./workflow/scripts/install-env.sh          # into $HOME/software/aligner-index-builder-env
source workflow/scripts/activate-env.sh    # activate in the current shell
```

`install-env.sh` downloads `aligner-index-builder-<version>-env.tar.gz` (+
`SHA256SUMS`) from the latest GitHub Release and unpacks it. Options:

- `-o PREFIX` -- install elsewhere (e.g. shared storage for a cluster)
- `-r vX.Y.Z` -- pin a specific release instead of latest
- `-f` -- overwrite even if `PREFIX` holds unrelated files

On SLURM compute nodes you can skip activation and just use the binaries:
`export PATH="$PREFIX/bin:$PATH"`.

### Option B -- conda (any platform)

```bash
conda env create -f workflow/environment.yaml
conda activate aligner-index-builder
```

## Usage

```bash
python workflow/scripts/aligner-index-builder \
  --fasta genome.fa --gtf genes.gtf \
  -o results/GRCh38 \
  --tools all \
  --cores 16
```

`--fasta` and `--gtf` may be local paths (plain, `.gz` or `.bz2` -- the
compression is sniffed, not assumed) or `http(s)://`, `ftp://`, `file://`
URLs. The CLI:

1. resolves local paths to absolute paths,
2. writes the run config to `<outdir>/config/<genome_name>.config.yaml`,
3. runs snakemake on the bundled `workflow/Snakefile`,
4. on success writes a manifest to `<outdir>/config/<genome_name>.yaml`.

A run always produces the FASTA (normalized + decompressed) with `.fai` and
`.2bit`, so you can hand those on to downstream pipelines.

### CLI options

| Option | Description |
| --- | --- |
| `--fasta PATH` | reference FASTA (required); local path or URL |
| `--gtf PATH` | gene annotation GTF (local path or URL); required for `star`, `hisat2`, `salmon`, which are otherwise skipped |
| `-o, --outdir DIR` | output directory, created if missing (required) |
| `--genome-name NAME` | name for the generated config + manifest (default: outdir basename) |
| `--tools TOOL...` | space- or comma-separated tools, or `all` (default) |
| `--sjdb-overhang N` | STAR splice-junction overhang, usually read length − 1 (default 100) |
| `--star-extra FLAGS` | extra `STAR --runMode genomeGenerate` flags (e.g. small-genome flags, see [Caveats](#caveats)) |
| `--salmon-kmer N` | salmon k-mer length (default 31) |
| `--cores N` | snakemake `--cores` (default 1) |
| `--slurm` | submit jobs to SLURM via the bundled profile `workflow/profiles/slurm` (see [HPC / SLURM](#hpc--slurm)) |
| `--dry-run` | write the config and print the snakemake plan, run nothing |
| `--unlock` | clear a stale snakemake lock left by an interrupted run (reuses `<outdir>/config/<name>.config.yaml`; only `-o`/`--genome-name` needed) |
| `--keep-temp` | keep intermediates snakemake would delete |
| `--snakemake-options ARGS` | extra snakemake arguments, passed through verbatim (quoted) |
| `-v, --verbose` | print snakemake's shell commands |
| `--version` | print version |

## HPC / SLURM

The repo ships a bundled [Snakemake workflow profile](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html)
for SLURM clusters (`workflow/profiles/slurm/config.yaml`). Pass `--slurm`
to the CLI and every job is submitted with `sbatch` instead of running
locally:

```bash
python workflow/scripts/aligner-index-builder \
  --fasta genome.fa --gtf genes.gtf \
  -o results/GRCh38 --tools all --slurm
```

The same profile works for direct snakemake runs (from the repo root, passing
`--directory <outdir>` so `.snakemake` stays in the output dir):

```bash
snakemake --configfile results/GRCh38/config/GRCh38.config.yaml \
  --workflow-profile workflow/profiles/slurm --directory results/GRCh38 \
  --cores 64
snakemake --configfile results/GRCh38/config/GRCh38.config.yaml \
  --workflow-profile workflow/profiles/slurm --directory results/GRCh38 \
  -n   # dry-run, nothing submitted
```

Metadata and locking live in the output directory, not the repo:

- Every run keeps its `.snakemake` (job metadata + the run lock) inside
  `<outdir>/.snakemake`, so the repo clone stays pristine and different
  output directories never contend for the same lock.
- If a run is interrupted — an `scancel`, a compute-node timeout, or a
  dropped ssh session on a `--slurm` run — snakemake leaves its lock behind
  and the next run fails with "directory is locked". Clear it with the CLI's
  `--unlock` (reuses the run config already written into the outdir, so only
  `-o` is needed):

  ```bash
  aligner-index-builder --unlock -o results/GRCh38
  ```

How jobs are sized:

- Each rule's `threads`, `mem_mb` and `runtime` come from
  `workflow/default-config/resources.yaml`. For a real genome the defaults
  are tuned to ~80% RAM efficiency against a measured GRCm38 run (index
  builds use little CPU, so all 9 index rules request 2 threads); add a
  `resources:` block to the generated config (deep-merged over the defaults)
  to override e.g. STAR's 40 GB, or pass
  `--snakemake-options '--default-resources mem_mb=... runtime=...'`.
- The profile's own `default-resources` only cover rules that declare none.
- `hisat2` is the big one: the splice-aware build (`hisat2-build
  --ss/--exon`, embedding the annotated splice sites and exons into the
  index) needs ~190 GB RAM for human-sized genomes (measured ~146 GB on
  GRCm38), so it requests a 190 GiB node (`mem_mb: 194560`). If your cluster
  has no big-mem node, build the plain index instead (~8-16 GB): in
  `workflow/rules/index.smk`, remove the
  `--ss ... --exon ...` flags and the `splicesites:`/`exons:` inputs from the
  `hisat2` rule, and lower its memory in the generated config:

  ```yaml
  resources:
    hisat2:
      mem_mb: 16000
  ```

  HISAT2 then finds splice junctions de novo at alignment time. Note that
  this edits the workflow source, so a later `git pull` reverts it.
- `bwa_mem2`/`bwameth2` can rival `hisat2`: bwa-mem2's index build needs
  ~50-90 GB for a human-sized genome and scales with the reference size, so
  the defaults request 80 GiB (`bwa_mem2`, measured ~61 GB on GRCm38) and
  160 GiB (`bwameth2`, measured ~122 GB). An
  undersized allocation is OOM-killed mid-build and `bwameth.py` surfaces it
  as `return code was:-9` inside an otherwise-empty rule log — raise
  `resources:` (or drop `threads` for the `bwa_mem2`/`bwameth2` rules, which
  lowers peak RAM) rather than reading it as a tool error. Smaller genomes
  can lower these in the generated config.

Cluster specifics to edit in `workflow/profiles/slurm/config.yaml`:

- `slurm_account` (`icmm_dm`) and `qos` (`normal`) are the ICMM group values
  — replace them if you are on another cluster or in another group.
- `slurm_partition` is intentionally unset so `sbatch` uses the cluster's
  default partition; set it only if yours has no default.
- `latency-wait: 60` gives NFS time to show up after STAR/Bismark finish
  (they write `directory()` outputs); bump it if you ever see spurious
  "missing files after X seconds" failures.

Compute nodes need the tools: the workflow runs without `--sdm conda` (each
job inherits the submit environment), so the pre-built environment must be
visible there — prepend its `bin` to `PATH` on the nodes, e.g. with
`export PATH="$PREFIX/bin:$PATH"` (see [Install](#install)).

`--slurm` also works with `--dry-run` (validates the profile, submits
nothing). `sbatch` must be available wherever snakemake runs.

## What gets built

| Tool | Index dir | Marker file | Notes |
| --- | --- | --- | --- |
| STAR | `index/star/` | `SAindex` | needs GTF (`--sjdbGTFfile`, overhang from `--sjdb-overhang`) |
| HISAT2 | `index/hisat2/` | `genome.6.ht2` | needs GTF (splice sites + exons extracted with `hisat2_extract_splice_sites.py`/`_exons.py`) |
| bowtie2 | `index/bowtie2/` | `genome.rev.2.bt2` | FASTA only |
| bwa | `index/bwa/` | `genome.fa.sa` | FASTA only |
| bwa-mem2 | `index/bwa-mem2/` | `genome.fa.bwt.2bit.64` | FASTA only |
| bwameth | `index/bwameth/` | `genome.fa.bwameth.c2t.sa` | FASTA only; `bwameth.py index` |
| bwameth2 | `index/bwameth2/` | `genome.fa.bwameth.c2t.bwt.2bit.64` | FASTA only; `bwameth.py index-mem2` |
| Bismark | `index/bismark/` | `Bisulfite_Genome/CT_conversion/BS_CT.1.bt2` | FASTA only; bowtie2-based (`bismark_genome_preparation --bowtie2`). Point `bismark --genome` at `index/bismark/` (the parent, not the `Bisulfite_Genome/` subdir) |
| salmon | `index/salmon/` | `seq.bin` | needs GTF; decoy-aware (transcripts via `gffread` + genome as decoys) |
| — | `genome_fasta/` | `genome.fa`, `genome.fa.fai`, `genome.2bit` | always; `samtools faidx` + `faToTwoBit` |
| — | `annotation/` | `genes.gtf`, `transcripts.fa` | when a GTF is given |
| — | `qc/` | `multiqc_report.html`, `multiqc_report_data/` | always; MultiQC report (see [MultiQC report](#multiqc-report)) |

Tool version pins live in `workflow/default-config/versions.yaml` (also
mirrored into `workflow/environment.yaml`). Note: STAR 2.7.11b is the final
STAR release and pins `htslib < 1.23`, so `samtools` is pinned to `1.22.1` to
stay in the same environment.

## Output layout

```
results/GRCh38/
├── config/
│   ├── GRCh38.config.yaml       # run config (written by the CLI)
│   └── GRCh38.yaml              # success manifest (written by the CLI)
├── genome_fasta/
│   ├── genome.fa               # normalized, uncompressed reference
│   ├── genome.fa.fai
│   └── genome.2bit
├── annotation/
│   ├── genes.gtf               # normalized annotation
│   └── transcripts.fa          # gffread transcriptome
├── index/
│   ├── star/
│   ├── hisat2/
│   ├── bowtie2/
│   ├── bwa/
│   ├── bwa-mem2/
│   ├── bwameth/
│   ├── bwameth2/
│   ├── bismark/Bisulfite_Genome/
│   └── salmon/
├── versions/
│   └── GRCh38_mqc_versions.yml     # tool versions for the report
├── qc/
│   ├── multiqc_report.html         # the MultiQC report
│   └── multiqc_report_data/        # machine-readable report data
└── pipeline_info/
    ├── logs/                       # per-rule logs + config_resolution.log
    ├── benchmarks/                 # per-rule timing/peak-memory tables
    ├── chrom_consistency.json      # FASTA-vs-GTF contig check result
    ├── resources_used_mqc.json     # report: "Resources used" table
    ├── index_sizes_mqc.json        # report: "Index sizes" table
    ├── annotation_summary_mqc.json # report: "Annotation summary" table (with GTF)
    ├── config_used_mqc.json        # report: "Configuration used" table
    └── benchmark_summary_mqc.json  # report: "Resource usage" table
```

## MultiQC report

The pipeline has no per-sample data (no reads, no alignments), so the report in
`qc/multiqc_report.html` is custom content aggregated by MultiQC from the files
under `versions/` and `pipeline_info/`:

- **Resources used** -- reference FASTA path + SHA-256, contig count / total
  length / N50 / GC%, GTF path + SHA-256, the FASTA-vs-GTF contig-consistency
  result, and the index sets built (`resources_used_mqc.json`).
- **Index sizes** -- disk footprint of each built index directory
  (`index_sizes_mqc.json`).
- **Annotation summary** -- gene/transcript counts from the GTF (+ the gffread
  transcriptome when salmon is built). Omitted when no GTF is given
  (`annotation_summary_mqc.json`).
- **Configuration used** -- the resolved run config (`config_used_mqc.json`).
- **Resource usage** -- per-rule job count, wall time, and CPU/RAM efficiency
  (allocated vs actually used, from the Snakemake benchmark tables)
  (`benchmark_summary_mqc.json`). Handy for sizing your cluster before a full run.
- **Software Versions** -- the pinned tool versions from
  `versions.yaml` (`versions/<genome>_mqc_versions.yml`, auto-discovered).

The report title, comment and section order live in
`workflow/default-config/multiqc_config.yaml`. The report is cheap (seconds),
so it is always produced.

## Direct snakemake use

For debugging or scripting, generate the config with `--dry-run`, or write one
by hand (see `config/config.example.yaml`) and run:

```bash
snakemake --configfile results/GRCh38/config/GRCh38.config.yaml --cores N
```

CLI-generated configs carry a `repo_root` key, which anchors the
`workflow/default-config/*.yaml` includes regardless of the working
directory (the CLI runs snakemake with `--directory <outdir>`). For a
hand-written config, either add `repo_root` too or run snakemake from the
repo root.

The config chain loaded by `workflow/Snakefile` is:
`workflow/default-config/{versions,resources}.yaml` (built-in) &rarr; your run
config &rarr; optional `resources:` override block inside the run config.

## Config reference

| Key | Type | Default | Meaning |
| --- | --- | --- | --- |
| `genome_name` | string | outdir basename | name for config/manifest |
| `outdir` | string | — | output directory (required) |
| `repo_root` | string | cwd | path to the `aligner-index-builder` repo (written by the CLI; anchors the `default-config` includes) |
| `fasta` | string | — | FASTA path or URL (required) |
| `gtf` | string/null | `null` | GTF path or URL; enables `star`/`hisat2`/`salmon` + the chromosome-consistency check |
| `tools` | array or `all` | — | which index sets to build (required) |
| `sjdb_overhang` | int | `100` | STAR `--sjdbOverhang` |
| `star_extra` | string | `""` | extra STAR genomeGenerate flags |
| `salmon_kmer` | int | `31` | salmon `-k` |
| `resources` | map | — | per-rule `{threads, mem_mb, runtime}` overrides (deep-merged over `default-config/resources.yaml`) |

## Architecture

```
aligner-index-builder (CLI)
   └─ writes <outdir>/config/<name>.config.yaml
       └─ snakemake -s workflow/Snakefile --configfile <config>
           ├─ default-config/{versions,resources,multiqc_config}.yaml
           │                      (built-in defaults)
           ├─ rules/common.smk  (config validation, paths, per-tool envs)
           ├─ rules/index.smk   (one rule per index)
           │   ├─ fetch_reference.py (download/copy + gzip/bz2 sniff-decompress)
           │   └─ check_chroms.py   (FASTA-vs-GTF contig consistency gate)
           └─ rules/qc.smk      (MultiQC report; see "MultiQC report")
               ├─ resources_used.py    (reference/annotation/index-size tables)
               └─ benchmark_summary.py (per-rule resource-usage table)
```

```mermaid
flowchart LR
    fasta[(FASTA)] --> prep[prepare_genome]
    gtf[(GTF)] --> prepgtf[prepare_gtf]
    prep --> fai[samtools faidx]
    prep --> twoBit[faToTwoBit]
    prep --> star[STAR index]
    prep --> hisat2[HISAT2 index]
    prep --> bowtie2[bowtie2 index]
    prep --> bwa[bwa index]
    prep --> bwa2[bwa-mem2 index]
    prep --> bwameth[bwameth index]
    prep --> bwameth2[bwameth2 index]
    prep --> bismark[Bismark index]
    prep --> salmon[salmon index]
    prepgtf --> check[check_chrom_consistency]
    prepgtf --> star
    prepgtf --> hisat2
    prepgtf --> gffread[gffread transcripts]
    gffread --> salmon
    check --> qc[multiqc]
    star --> qc
    hisat2 --> qc
    bowtie2 --> qc
    bwa --> qc
    bwa2 --> qc
    bwameth --> qc
    bwameth2 --> qc
    bismark --> qc
    salmon --> qc
```

Per-tool conda env files are generated into `workflow/envs/generated/` at
parse time from `versions.yaml`; they're used only with
`snakemake --sdm conda`. The pre-built env already ships the same tools, so
plain `snakemake` runs use the activated environment directly.

## Versioning & releases

- `VERSION` holds the project version; `versions.yaml` holds tool versions.
- Releases are tagged `vX.Y.Z`. Pushing a tag triggers
  `.github/workflows/release-env.yml`, which conda-packs the environment,
  splits it into `<2 GiB` parts if needed, writes `SHA256SUMS`, and creates
  the release. `install-env.sh` resolves the latest release (API-free: the
  releases Atom feed + the shipped `SHA256SUMS`; `AIB_INSTALL_ORIGIN`
  overrides the origin host for mirrors/tests).
- CI (`.github/workflows/ci.yml`) builds every index set from the tiny
  committed fixture (`.tests/`) and runs the negative edge-case guards.

## Caveats

- **STAR on small genomes**: STAR aborts `genomeGenerate` on very small
  references unless told otherwise. For the ~120 kb CI fixture the flag is
  `--star-extra "--genomeSAindexNbases 7 --genomeChrBinNbits 16"` (printed by
  `.tests/generate_test_data.py`); real chromosomes are unaffected.
- **STAR exit code**: STAR has a long-standing benign exit-time crash (a
  segfault *after* all index files are written). The `star` rule treats a
  non-zero exit as success if the core files actually landed on disk.
- **macOS**: the pre-built env is Linux x86_64 only. Build it locally with
  conda (Option B) or run on CI/linux.
- **Reproducibility**: tool builds are deterministic for a given reference;
  outputs land in `outdir` and are gitignored, so the repo stays a clean
  project and `git pull` never conflicts with generated data.
