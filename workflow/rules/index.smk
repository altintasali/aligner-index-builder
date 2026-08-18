# =============================================================================
# Index-building rules. Each tool gets its own rule and writes into its own
# directory under <outdir>/index/<tool>/ (index/star/, index/hisat2/,
# index/bowtie2/, index/bwa/, index/bwa-mem2/, index/bwameth/,
# index/bwameth2/, index/bismark/, index/salmon/) plus the shared
# genome_fasta/ and annotation/ side-products.
#
# Rules whose tool needs the GTF (STAR, HISAT2, salmon) are defined only when
# a GTF was given -- common.smk drops those tools from TOOLS otherwise.
# =============================================================================

# -----------------------------------------------------------------------------
# Reference preparation
# -----------------------------------------------------------------------------
rule prepare_genome:
    # Copies a local FASTA (decompressing it if gzipped/bzipped) or downloads
    # a URL into outdir/genome_fasta/genome_raw.fa.  No input file is declared:
    # the source may be a URL, which snakemake can't treat as an input file.
    # The raw FASTA is always produced; filter_genome (the next step) either
    # copies it verbatim or subsets it based on --keep-chroms / --remove-
    # nonstandard / --exclude-chroms before any index is built.
    output:
        GENOME_RAW_FASTA,
    params:
        src=FASTA,
    threads: get_resources("prepare_genome")["threads"]
    resources:
        mem_mb=get_resources("prepare_genome")["mem_mb"],
        runtime=get_resources("prepare_genome")["runtime"],
    benchmark:
        os.path.join(OUTDIR, "pipeline_info", "benchmarks", "prepare_genome.txt"),
    log:
        os.path.join(OUTDIR, "pipeline_info", "logs", "prepare_genome.log"),
    conda:
        PYTHON_ENV
    script:
        "../scripts/fetch_reference.py"


rule filter_genome:
    # Subsets the raw FASTA based on chromosome name rules.  When no filtering
    # is requested the script copies the file verbatim so the DAG is the same
    # either way -- downstream rules always read GENOME_FASTA.
    input:
        GENOME_RAW_FASTA,
    output:
        GENOME_FASTA,
    params:
        remove_nonstandard=REMOVE_NONSTANDARD,
        keep_chroms=KEEP_CHROMS,
        exclude_chroms=EXCLUDE_CHROMS,
    threads: get_resources("filter_genome")["threads"]
    resources:
        mem_mb=get_resources("filter_genome")["mem_mb"],
        runtime=get_resources("filter_genome")["runtime"],
    benchmark:
        os.path.join(OUTDIR, "pipeline_info", "benchmarks", "filter_genome.txt"),
    log:
        os.path.join(OUTDIR, "pipeline_info", "logs", "filter_genome.log"),
    conda:
        PYTHON_ENV
    script:
        "../scripts/filter_genome.py"


rule fasta_fai:
    input:
        GENOME_FASTA,
    output:
        GENOME_FAI,
    threads: get_resources("fasta_fai")["threads"]
    resources:
        mem_mb=get_resources("fasta_fai")["mem_mb"],
        runtime=get_resources("fasta_fai")["runtime"],
    benchmark:
        os.path.join(OUTDIR, "pipeline_info", "benchmarks", "fasta_fai.txt"),
    log:
        os.path.join(OUTDIR, "pipeline_info", "logs", "fasta_fai.log"),
    conda:
        SAMTOOLS_ENV
    shell:
        "samtools faidx {input}"


rule fasta_2bit:
    input:
        GENOME_FASTA,
    output:
        GENOME_2BIT,
    threads: get_resources("fasta_2bit")["threads"]
    resources:
        mem_mb=get_resources("fasta_2bit")["mem_mb"],
        runtime=get_resources("fasta_2bit")["runtime"],
    benchmark:
        os.path.join(OUTDIR, "pipeline_info", "benchmarks", "fasta_2bit.txt"),
    log:
        os.path.join(OUTDIR, "pipeline_info", "logs", "fasta_2bit.log"),
    conda:
        UCSC_ENV
    shell:
        "faToTwoBit {input} {output}"


# -----------------------------------------------------------------------------
# STAR (needs the GTF for splice junctions)
# -----------------------------------------------------------------------------
rule star:
    # STAR has a long-standing, benign crash-on-exit bug (a segfault *after*
    # all output is written) -- see alexdobin/STAR issues; confirmed by the
    # STAR author not to affect results. So on a non-zero exit we check
    # whether the core index files actually landed on disk before treating it
    # as a real failure. The CLI's --star-extra is where you'd set the
    # small-genome flags (--genomeSAindexNbases/--genomeChrBinNbits).
    input:
        fasta=GENOME_FASTA,
        gtf=GENES_GTF if GENES_GTF else [],
    output:
        directory(os.path.join(INDEX_DIR, "star")),
    params:
        gtf_arg=(f"--sjdbGTFfile {GENES_GTF} " if GENES_GTF else ""),
        extra=config.get("star_extra", "--sjdbOverhang 100"),
    threads: get_resources("star")["threads"]
    resources:
        mem_mb=get_resources("star")["mem_mb"],
        runtime=get_resources("star")["runtime"],
    benchmark:
        os.path.join(OUTDIR, "pipeline_info", "benchmarks", "star.txt"),
    log:
        os.path.join(OUTDIR, "pipeline_info", "logs", "star.log"),
    conda:
        STAR_ENV
    shell:
        "mkdir -p {output} && "
        "(cd {output} && "
        "STAR --runMode genomeGenerate "
        "--genomeDir {output} "
        "--genomeFastaFiles {input.fasta} "
        "{params.gtf_arg}"
        "--runThreadN {threads} "
        "{params.extra}) "
        "> {log} 2>&1 "
        "|| (echo 'STAR exited non-zero; checking whether the index was "
        "actually written successfully anyway (benign STAR exit-time crash)' "
        ">> {log}; "
        "test -s {output}/SA && test -s {output}/SAindex && "
        "test -s {output}/Genome)"


# -----------------------------------------------------------------------------
# HISAT2 (uses GTF-derived splice sites and exons)
# -----------------------------------------------------------------------------
if GENES_GTF:

    rule hisat2_splicesites:
        input:
            GENES_GTF,
        output:
            temp(os.path.join(INDEX_DIR, "hisat2", "splice_sites.tsv")),
        threads: get_resources("hisat2_splicesites")["threads"]
        resources:
            mem_mb=get_resources("hisat2_splicesites")["mem_mb"],
            runtime=get_resources("hisat2_splicesites")["runtime"],
        benchmark:
            os.path.join(OUTDIR, "pipeline_info", "benchmarks", "hisat2_splicesites.txt"),
        log:
            os.path.join(OUTDIR, "pipeline_info", "logs", "hisat2_splicesites.log"),
        conda:
            HISAT2_ENV
        shell:
            "mkdir -p {INDEX_DIR}/hisat2 && "
            "hisat2_extract_splice_sites.py {input} > {output}"

    rule hisat2_exons:
        input:
            GENES_GTF,
        output:
            temp(os.path.join(INDEX_DIR, "hisat2", "exons.tsv")),
        threads: get_resources("hisat2_exons")["threads"]
        resources:
            mem_mb=get_resources("hisat2_exons")["mem_mb"],
            runtime=get_resources("hisat2_exons")["runtime"],
        benchmark:
            os.path.join(OUTDIR, "pipeline_info", "benchmarks", "hisat2_exons.txt"),
        log:
            os.path.join(OUTDIR, "pipeline_info", "logs", "hisat2_exons.log"),
        conda:
            HISAT2_ENV
        shell:
            "mkdir -p {INDEX_DIR}/hisat2 && "
            "hisat2_extract_exons.py {input} > {output}"

    rule hisat2:
        input:
            fasta=GENOME_FASTA,
            splicesites=os.path.join(INDEX_DIR, "hisat2", "splice_sites.tsv"),
            exons=os.path.join(INDEX_DIR, "hisat2", "exons.tsv"),
        output:
            os.path.join(INDEX_DIR, "hisat2", "genome.6.ht2"),
        params:
            extra=config.get("hisat2_extra", ""),
        threads: get_resources("hisat2")["threads"]
        resources:
            mem_mb=get_resources("hisat2")["mem_mb"],
            runtime=get_resources("hisat2")["runtime"],
        benchmark:
            os.path.join(OUTDIR, "pipeline_info", "benchmarks", "hisat2.txt"),
        log:
            os.path.join(OUTDIR, "pipeline_info", "logs", "hisat2.log"),
        conda:
            HISAT2_ENV
        shell:
            "hisat2-build -q -p {threads} "
            "{params.extra} "
            "--ss {input.splicesites} "
            "--exon {input.exons} "
            "{input.fasta} {INDEX_DIR}/hisat2/genome"


# -----------------------------------------------------------------------------
# bowtie2 (FASTA only)
# -----------------------------------------------------------------------------
rule bowtie2:
    input:
        GENOME_FASTA,
    output:
        os.path.join(INDEX_DIR, "bowtie2", "genome.rev.2.bt2"),
    params:
        extra=config.get("bowtie2_extra", ""),
    threads: get_resources("bowtie2")["threads"]
    resources:
        mem_mb=get_resources("bowtie2")["mem_mb"],
        runtime=get_resources("bowtie2")["runtime"],
    benchmark:
        os.path.join(OUTDIR, "pipeline_info", "benchmarks", "bowtie2.txt"),
    log:
        os.path.join(OUTDIR, "pipeline_info", "logs", "bowtie2.log"),
    conda:
        BOWTIE2_ENV
    shell:
        "mkdir -p {INDEX_DIR}/bowtie2 && "
        "bowtie2-build --threads {threads} {params.extra} "
        "{input} {INDEX_DIR}/bowtie2/genome"


# -----------------------------------------------------------------------------
# bwa / bwa-mem2 (FASTA only). The fasta is symlinked into the index dir so
# each index is self-contained while avoiding a full copy of a multi-GB file.
# -----------------------------------------------------------------------------
rule bwa:
    input:
        GENOME_FASTA,
    output:
        os.path.join(INDEX_DIR, "bwa", "genome.fa.sa"),
    params:
        extra=config.get("bwa_extra", ""),
    threads: get_resources("bwa")["threads"]
    resources:
        mem_mb=get_resources("bwa")["mem_mb"],
        runtime=get_resources("bwa")["runtime"],
    benchmark:
        os.path.join(OUTDIR, "pipeline_info", "benchmarks", "bwa.txt"),
    log:
        os.path.join(OUTDIR, "pipeline_info", "logs", "bwa.log"),
    conda:
        BWA_ENV
    shell:
        "mkdir -p {INDEX_DIR}/bwa && "
        "ln -sf {input} {INDEX_DIR}/bwa/genome.fa && "
        "bwa index {params.extra} {INDEX_DIR}/bwa/genome.fa"


rule bwa_mem2:
    input:
        GENOME_FASTA,
    output:
        os.path.join(INDEX_DIR, "bwa-mem2", "genome.fa.bwt.2bit.64"),
    params:
        extra=config.get("bwa_mem2_extra", ""),
    threads: get_resources("bwa_mem2")["threads"]
    resources:
        mem_mb=get_resources("bwa_mem2")["mem_mb"],
        runtime=get_resources("bwa_mem2")["runtime"],
    benchmark:
        os.path.join(OUTDIR, "pipeline_info", "benchmarks", "bwa_mem2.txt"),
    log:
        os.path.join(OUTDIR, "pipeline_info", "logs", "bwa_mem2.log"),
    conda:
        BWA_MEM2_ENV
    shell:
        "mkdir -p {INDEX_DIR}/bwa-mem2 && "
        "ln -sf {input} {INDEX_DIR}/bwa-mem2/genome.fa && "
        "bwa-mem2 index {params.extra} {INDEX_DIR}/bwa-mem2/genome.fa"


# -----------------------------------------------------------------------------
# bwameth / bwameth2 (bisulfite-aware, FASTA only)
# -----------------------------------------------------------------------------
rule bwameth:
    input:
        GENOME_FASTA,
    output:
        os.path.join(INDEX_DIR, "bwameth", "genome.fa.bwameth.c2t.sa"),
    params:
        extra=config.get("bwameth_extra", ""),
    threads: get_resources("bwameth")["threads"]
    resources:
        mem_mb=get_resources("bwameth")["mem_mb"],
        runtime=get_resources("bwameth")["runtime"],
    benchmark:
        os.path.join(OUTDIR, "pipeline_info", "benchmarks", "bwameth.txt"),
    log:
        os.path.join(OUTDIR, "pipeline_info", "logs", "bwameth.log"),
    conda:
        BWAMETH_ENV
    shell:
        "mkdir -p {INDEX_DIR}/bwameth && "
        "ln -sf {input} {INDEX_DIR}/bwameth/genome.fa && "
        "bwameth.py index {params.extra} {INDEX_DIR}/bwameth/genome.fa"


rule bwameth2:
    input:
        GENOME_FASTA,
    output:
        os.path.join(INDEX_DIR, "bwameth2", "genome.fa.bwameth.c2t.bwt.2bit.64"),
    params:
        extra=config.get("bwameth2_extra", ""),
    threads: get_resources("bwameth2")["threads"]
    resources:
        mem_mb=get_resources("bwameth2")["mem_mb"],
        runtime=get_resources("bwameth2")["runtime"],
    benchmark:
        os.path.join(OUTDIR, "pipeline_info", "benchmarks", "bwameth2.txt"),
    log:
        os.path.join(OUTDIR, "pipeline_info", "logs", "bwameth2.log"),
    conda:
        BWAMETH_ENV
    shell:
        "mkdir -p {INDEX_DIR}/bwameth2 && "
        "ln -sf {input} {INDEX_DIR}/bwameth2/genome.fa && "
        "bwameth.py index-mem2 {params.extra} {INDEX_DIR}/bwameth2/genome.fa"


# -----------------------------------------------------------------------------
# Bismark (bisulfite, FASTA only). bismark_genome_preparation builds the
# bowtie2-based Bisulfite_Genome/ (CT + GA conversion) inside the index dir.
# -----------------------------------------------------------------------------
rule bismark:
    input:
        GENOME_FASTA,
    output:
        directory(os.path.join(INDEX_DIR, "bismark", "Bisulfite_Genome")),
    params:
        extra=config.get("bismark_extra", ""),
    threads: get_resources("bismark")["threads"]
    resources:
        mem_mb=get_resources("bismark")["mem_mb"],
        runtime=get_resources("bismark")["runtime"],
    benchmark:
        os.path.join(OUTDIR, "pipeline_info", "benchmarks", "bismark.txt"),
    log:
        os.path.join(OUTDIR, "pipeline_info", "logs", "bismark.log"),
    conda:
        BISMARK_ENV
    shell:
        "mkdir -p {INDEX_DIR}/bismark && "
        "ln -sf {input} {INDEX_DIR}/bismark/genome.fa && "
        "bismark_genome_preparation --bowtie2 "
        "{params.extra} "
        "--path_to_aligner \"$(dirname \"$(command -v bowtie2)\")\" "
        "{INDEX_DIR}/bismark"


# -----------------------------------------------------------------------------
# GTF-derived annotation + chromosome-consistency check (needs the GTF)
# -----------------------------------------------------------------------------
if GENES_GTF:

    rule prepare_gtf:
        output:
            GENES_GTF,
        params:
            src=GTF,
        threads: get_resources("prepare_gtf")["threads"]
        resources:
            mem_mb=get_resources("prepare_gtf")["mem_mb"],
            runtime=get_resources("prepare_gtf")["runtime"],
        benchmark:
            os.path.join(OUTDIR, "pipeline_info", "benchmarks", "prepare_gtf.txt"),
        log:
            os.path.join(OUTDIR, "pipeline_info", "logs", "prepare_gtf.log"),
        conda:
            PYTHON_ENV
        script:
            "../scripts/fetch_reference.py"

    rule check_chrom_consistency:
        # Fails fast with a clear message if the FASTA and GTF don't share at
        # least one chromosome/contig name -- otherwise STAR/HISAT2/salmon
        # would build indices against a mismatched reference. Also writes the
        # shared-contig counts consumed by the report's resources_used rule.
        input:
            fasta=GENOME_FASTA,
            gtf=GENES_GTF,
        output:
            ok=touch(CHROM_CHECK_OK),
            stats=CHROM_STATS,
        threads: get_resources("check_chrom_consistency")["threads"]
        resources:
            mem_mb=get_resources("check_chrom_consistency")["mem_mb"],
            runtime=get_resources("check_chrom_consistency")["runtime"],
        benchmark:
            os.path.join(OUTDIR, "pipeline_info", "benchmarks", "chrom_consistency.txt"),
        log:
            os.path.join(OUTDIR, "pipeline_info", "logs", "chrom_consistency.log"),
        conda:
            PYTHON_ENV
        script:
            "../scripts/check_chroms.py"


# -----------------------------------------------------------------------------
# Salmon (decoy-aware: transcriptome from gffread + genome appended as decoys)
# -----------------------------------------------------------------------------
if GENES_GTF:

    rule transcriptome_fasta:
        input:
            fasta=GENOME_FASTA,
            gtf=GENES_GTF,
        output:
            TRANSCRIPTS_FA,
        threads: get_resources("transcriptome_fasta")["threads"]
        resources:
            mem_mb=get_resources("transcriptome_fasta")["mem_mb"],
            runtime=get_resources("transcriptome_fasta")["runtime"],
        benchmark:
            os.path.join(OUTDIR, "pipeline_info", "benchmarks", "transcriptome_fasta.txt"),
        log:
            os.path.join(OUTDIR, "pipeline_info", "logs", "transcriptome_fasta.log"),
        conda:
            GFFREAD_ENV
        shell:
            "gffread -w {output} -g {input.fasta} {input.gtf}"

    rule salmon:
        input:
            fasta=GENOME_FASTA,
            transcripts=TRANSCRIPTS_FA,
        output:
            seq_bin=os.path.join(INDEX_DIR, "salmon", "seq.bin"),
            seq_fa=temp(os.path.join(INDEX_DIR, "salmon", "seq.fa")),
            decoys=os.path.join(INDEX_DIR, "salmon", "decoys.txt"),
        params:
            extra=config.get("salmon_extra", "-k 31"),
        threads: get_resources("salmon")["threads"]
        resources:
            mem_mb=get_resources("salmon")["mem_mb"],
            runtime=get_resources("salmon")["runtime"],
        benchmark:
            os.path.join(OUTDIR, "pipeline_info", "benchmarks", "salmon.txt"),
        log:
            os.path.join(OUTDIR, "pipeline_info", "logs", "salmon.log"),
        conda:
            SALMON_ENV
        # decoys.txt must hold only the first whitespace-delimited field of each
        # genome header and seq.fa must concatenate transcripts BEFORE the genome,
        # otherwise Salmon misidentifies chromosome IDs as transcript IDs at the
        # decoy boundary (salmon_desired_abs). Keep both invariants intact.
        shell:
            "mkdir -p {INDEX_DIR}/salmon && "
            "grep '^>' {input.fasta} | cut -d' ' -f1 | tr -d '>' > {output.decoys} && "
            "cat {input.transcripts} {input.fasta} > {output.seq_fa} && "
            "salmon index -p {threads} -t {output.seq_fa} "
            "-d {output.decoys} -i {INDEX_DIR}/salmon "
            "{params.extra}"
