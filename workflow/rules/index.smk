# =============================================================================
# Index-building rules. Each tool gets its own rule and writes into its own
# directory under OUTDIR (star_index/, hisat2_index/, bowtie2_index/,
# bwa_index/, bwa-mem2_index/, bwa-meth{2}_index/, bismark_index/,
# salmon_index/) plus the shared genome_fasta/ and annotation/ side-products.
#
# Rules whose tool needs the GTF (STAR, HISAT2, salmon) are defined only when
# a GTF was given -- common.smk drops those tools from TOOLS otherwise.
# =============================================================================

# -----------------------------------------------------------------------------
# Reference preparation
# -----------------------------------------------------------------------------
rule prepare_genome:
    # Copies a local FASTA (decompressing it if gzipped/bzipped) or downloads
    # a URL into outdir/genome_fasta/genome.fa. No input file is declared:
    # the source may be a URL, which snakemake can't treat as an input file.
    output:
        GENOME_FASTA,
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
rule star_index:
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
        directory(os.path.join(OUTDIR, "star_index")),
    params:
        sjdb_overhang=config.get("sjdb_overhang", 100),
        gtf_arg=(f"--sjdbGTFfile {GENES_GTF} " if GENES_GTF else ""),
        extra=config.get("star_extra", ""),
    threads: get_resources("star_index")["threads"]
    resources:
        mem_mb=get_resources("star_index")["mem_mb"],
        runtime=get_resources("star_index")["runtime"],
    benchmark:
        os.path.join(OUTDIR, "pipeline_info", "benchmarks", "star_index.txt"),
    log:
        os.path.join(OUTDIR, "pipeline_info", "logs", "star_index.log"),
    conda:
        STAR_ENV
    shell:
        "mkdir -p {output} && "
        "(cd {output} && "
        "STAR --runMode genomeGenerate "
        "--genomeDir {output} "
        "--genomeFastaFiles {input.fasta} "
        "{params.gtf_arg}"
        "--sjdbOverhang {params.sjdb_overhang} "
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
            os.path.join(OUTDIR, "hisat2_index", "splice_sites.tsv"),
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
            "mkdir -p {OUTDIR}/hisat2_index && "
            "hisat2_extract_splice_sites.py {input} > {output}"

    rule hisat2_exons:
        input:
            GENES_GTF,
        output:
            os.path.join(OUTDIR, "hisat2_index", "exons.tsv"),
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
            "mkdir -p {OUTDIR}/hisat2_index && "
            "hisat2_extract_exons.py {input} > {output}"

    rule hisat2_index:
        input:
            fasta=GENOME_FASTA,
            splicesites=os.path.join(OUTDIR, "hisat2_index", "splice_sites.tsv"),
            exons=os.path.join(OUTDIR, "hisat2_index", "exons.tsv"),
        output:
            os.path.join(OUTDIR, "hisat2_index", "genome.6.ht2"),
        threads: get_resources("hisat2_index")["threads"]
        resources:
            mem_mb=get_resources("hisat2_index")["mem_mb"],
            runtime=get_resources("hisat2_index")["runtime"],
        benchmark:
            os.path.join(OUTDIR, "pipeline_info", "benchmarks", "hisat2_index.txt"),
        log:
            os.path.join(OUTDIR, "pipeline_info", "logs", "hisat2_index.log"),
        conda:
            HISAT2_ENV
        shell:
            "hisat2-build -q -p {threads} "
            "--ss {input.splicesites} "
            "--exon {input.exons} "
            "{input.fasta} {OUTDIR}/hisat2_index/genome"


# -----------------------------------------------------------------------------
# bowtie2 (FASTA only)
# -----------------------------------------------------------------------------
rule bowtie2_index:
    input:
        GENOME_FASTA,
    output:
        os.path.join(OUTDIR, "bowtie2_index", "genome.rev.2.bt2"),
    threads: get_resources("bowtie2_index")["threads"]
    resources:
        mem_mb=get_resources("bowtie2_index")["mem_mb"],
        runtime=get_resources("bowtie2_index")["runtime"],
    benchmark:
        os.path.join(OUTDIR, "pipeline_info", "benchmarks", "bowtie2_index.txt"),
    log:
        os.path.join(OUTDIR, "pipeline_info", "logs", "bowtie2_index.log"),
    conda:
        BOWTIE2_ENV
    shell:
        "mkdir -p {OUTDIR}/bowtie2_index && "
        "bowtie2-build --threads {threads} {input} {OUTDIR}/bowtie2_index/genome"


# -----------------------------------------------------------------------------
# bwa / bwa-mem2 (FASTA only). The fasta is symlinked into the index dir so
# each index is self-contained while avoiding a full copy of a multi-GB file.
# -----------------------------------------------------------------------------
rule bwa_index:
    input:
        GENOME_FASTA,
    output:
        os.path.join(OUTDIR, "bwa_index", "genome.fa.sa"),
    threads: get_resources("bwa_index")["threads"]
    resources:
        mem_mb=get_resources("bwa_index")["mem_mb"],
        runtime=get_resources("bwa_index")["runtime"],
    benchmark:
        os.path.join(OUTDIR, "pipeline_info", "benchmarks", "bwa_index.txt"),
    log:
        os.path.join(OUTDIR, "pipeline_info", "logs", "bwa_index.log"),
    conda:
        BWA_ENV
    shell:
        "mkdir -p {OUTDIR}/bwa_index && "
        "ln -sf {input} {OUTDIR}/bwa_index/genome.fa && "
        "bwa index {OUTDIR}/bwa_index/genome.fa"


rule bwa_mem2_index:
    input:
        GENOME_FASTA,
    output:
        os.path.join(OUTDIR, "bwa-mem2_index", "genome.fa.bwt.2bit.64"),
    threads: get_resources("bwa_mem2_index")["threads"]
    resources:
        mem_mb=get_resources("bwa_mem2_index")["mem_mb"],
        runtime=get_resources("bwa_mem2_index")["runtime"],
    benchmark:
        os.path.join(OUTDIR, "pipeline_info", "benchmarks", "bwa_mem2_index.txt"),
    log:
        os.path.join(OUTDIR, "pipeline_info", "logs", "bwa_mem2_index.log"),
    conda:
        BWA_MEM2_ENV
    shell:
        "mkdir -p {OUTDIR}/bwa-mem2_index && "
        "ln -sf {input} {OUTDIR}/bwa-mem2_index/genome.fa && "
        "bwa-mem2 index {OUTDIR}/bwa-mem2_index/genome.fa"


# -----------------------------------------------------------------------------
# bwameth / bwameth2 (bwa-meth; bisulfite-aware, FASTA only)
# -----------------------------------------------------------------------------
rule bwameth_index:
    input:
        GENOME_FASTA,
    output:
        os.path.join(OUTDIR, "bwa-meth_index", "genome.fa.bwameth.c2t.sa"),
    threads: get_resources("bwameth_index")["threads"]
    resources:
        mem_mb=get_resources("bwameth_index")["mem_mb"],
        runtime=get_resources("bwameth_index")["runtime"],
    benchmark:
        os.path.join(OUTDIR, "pipeline_info", "benchmarks", "bwameth_index.txt"),
    log:
        os.path.join(OUTDIR, "pipeline_info", "logs", "bwameth_index.log"),
    conda:
        BWAMETH_ENV
    shell:
        "mkdir -p {OUTDIR}/bwa-meth_index && "
        "ln -sf {input} {OUTDIR}/bwa-meth_index/genome.fa && "
        "bwameth.py index {OUTDIR}/bwa-meth_index/genome.fa"


rule bwameth2_index:
    input:
        GENOME_FASTA,
    output:
        os.path.join(OUTDIR, "bwa-meth2_index", "genome.fa.bwameth.c2t.bwt.2bit.64"),
    threads: get_resources("bwameth2_index")["threads"]
    resources:
        mem_mb=get_resources("bwameth2_index")["mem_mb"],
        runtime=get_resources("bwameth2_index")["runtime"],
    benchmark:
        os.path.join(OUTDIR, "pipeline_info", "benchmarks", "bwameth2_index.txt"),
    log:
        os.path.join(OUTDIR, "pipeline_info", "logs", "bwameth2_index.log"),
    conda:
        BWAMETH_ENV
    shell:
        "mkdir -p {OUTDIR}/bwa-meth2_index && "
        "ln -sf {input} {OUTDIR}/bwa-meth2_index/genome.fa && "
        "bwameth.py index-mem2 {OUTDIR}/bwa-meth2_index/genome.fa"


# -----------------------------------------------------------------------------
# Bismark (bisulfite, FASTA only). bismark_genome_preparation builds the
# bowtie2-based Bisulfite_Genome/ (CT + GA conversion) inside the index dir.
# -----------------------------------------------------------------------------
rule bismark_index:
    input:
        GENOME_FASTA,
    output:
        directory(os.path.join(OUTDIR, "bismark_index", "Bisulfite_Genome")),
    threads: get_resources("bismark_index")["threads"]
    resources:
        mem_mb=get_resources("bismark_index")["mem_mb"],
        runtime=get_resources("bismark_index")["runtime"],
    benchmark:
        os.path.join(OUTDIR, "pipeline_info", "benchmarks", "bismark_index.txt"),
    log:
        os.path.join(OUTDIR, "pipeline_info", "logs", "bismark_index.log"),
    conda:
        BISMARK_ENV
    shell:
        "mkdir -p {OUTDIR}/bismark_index && "
        "ln -sf {input} {OUTDIR}/bismark_index/genome.fa && "
        "bismark_genome_preparation --bowtie2 "
        "--path_to_aligner \"$(dirname \"$(command -v bowtie2)\")\" "
        "{OUTDIR}/bismark_index"


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

    rule salmon_index:
        input:
            fasta=GENOME_FASTA,
            transcripts=TRANSCRIPTS_FA,
        output:
            seq_bin=os.path.join(OUTDIR, "salmon_index", "seq.bin"),
            seq_fa=temp(os.path.join(OUTDIR, "salmon_index", "seq.fa")),
            decoys=os.path.join(OUTDIR, "salmon_index", "decoys.txt"),
        params:
            kmer=config.get("salmon_kmer", 31),
        threads: get_resources("salmon_index")["threads"]
        resources:
            mem_mb=get_resources("salmon_index")["mem_mb"],
            runtime=get_resources("salmon_index")["runtime"],
        benchmark:
            os.path.join(OUTDIR, "pipeline_info", "benchmarks", "salmon_index.txt"),
        log:
            os.path.join(OUTDIR, "pipeline_info", "logs", "salmon_index.log"),
        conda:
            SALMON_ENV
        shell:
            "mkdir -p {OUTDIR}/salmon_index && "
            "grep '^>' {input.fasta} | cut -d' ' -f1 | tr -d '>' > {output.decoys} && "
            "cat {input.transcripts} {input.fasta} > {output.seq_fa} && "
            "salmon index -p {threads} -t {output.seq_fa} "
            "-d {output.decoys} -i {OUTDIR}/salmon_index -k {params.kmer}"
