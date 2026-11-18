#!/usr/bin/env nextflow
// DAP-seq: trim -> align -> MAPQ/contig filter -> GEM peaks -> top-N peak sequences -> MEME-ChIP,
// then motif QC and gene annotation over all peaks, bigWig tracks + an IGV session, and a MultiQC report.
// Peaks are called per sample and, for replicate groups (samplesheet `group`), once more on the
// pooled replicates; everything after GEM runs on each of these peak sets.

include { PREP_GENOME; BOWTIE2_BUILD; TRIM_GALORE; BOWTIE2_ALIGN; SORT_FILTER_BAM; GEM; GEM_QC;
          TOP_PEAKS_FASTA; PEAK_QC; MEME_CHIP; MOTIF_QC; ANNOTATE_PEAKS; BIGWIG; IGV_SESSION;
          MULTIQC } from './modules.nf'

workflow {
    if (!params.input) error "Missing --input (samplesheet CSV)"
    if (!params.fasta) error "Missing --fasta (genome FASTA, plain or .gz)"
    if (!params.gtf)   error "Missing --gtf (gene models, plain or .gz, same contig names as --fasta)"

    def rows      = read_samplesheet(params.input)
    def keep      = keep_chromosomes()
    def motifs = qc_motifs_inputs()
    if (!motifs.consensus && !motifs.files && !params.qc_motifs_n_denovo)
        error "Motif QC has no motifs: set --qc_motifs or --qc_motifs_n_denovo > 0"

    // Each process gets the params it uses as an `opts` map, so -resume reruns it when they change
    // (see modules.nf header).
    // Reference (PREP_GENOME.out.fasta = the FASTA, decompressed if it was .gz)
    PREP_GENOME(file(params.fasta, checkIfExists: true), [keep: keep])
    def index
    if (params.bowtie2_index) {
        index = file(params.bowtie2_index, checkIfExists: true)
    } else {
        BOWTIE2_BUILD(PREP_GENOME.out.fasta)
        index = BOWTIE2_BUILD.out
    }

    // Reads -> filtered BAM (experiments and controls alike)
    def reads = channel.fromList(rows.collect { r ->
        def meta = [id: r.sample, single_end: !r.fastq_2, control: r.control ?: null]
        def fq   = [file(r.fastq_1, checkIfExists: true)]
        if (r.fastq_2) fq << file(r.fastq_2, checkIfExists: true)
        [meta, fq]
    })
    TRIM_GALORE(reads, [args: params.trim_args])
    BOWTIE2_ALIGN(TRIM_GALORE.out.reads, index, [args: params.bowtie2_args])
    SORT_FILTER_BAM(BOWTIE2_ALIGN.out.sam, [min_mapq: params.min_mapq, keep: keep])

    // Peak sets: each experiment, plus each replicate group. Gather each set's experiment and control
    // BAMs (lists, sorted by name so the command line, and the -resume key, is stable).
    def sets = peak_sets(rows)
    def size = sets.collectEntries { ps -> [ps.id, ps.expt.size() + ps.ctrl.size()] }
    def gem_in = channel.fromList(sets.collectMany { ps ->
            ps.expt.collect { s -> [s, ps.id, 'expt'] } + ps.ctrl.collect { s -> [s, ps.id, 'ctrl'] } })
        .combine(SORT_FILTER_BAM.out.filtered.map { meta, bam, bai -> [meta.id, bam, bai] }, by: 0)
        .map { _s, id, role, bam, bai -> [groupKey(id, size[id]), [role, bam, bai]] }
        .groupTuple()
        .map { id, files ->
            def sorted = files.sort { f -> f[1].name }
            def expt   = sorted.findAll { f -> f[0] == 'expt' }
            def ctrl   = sorted.findAll { f -> f[0] == 'ctrl' }
            [[id: id.toString()], expt*.getAt(1), expt*.getAt(2), ctrl*.getAt(1), ctrl*.getAt(2)]
        }

    GEM(gem_in, PREP_GENOME.out.sizes, PREP_GENOME.out.chr_dir,
        [jar: params.gem_jar, read_dist: params.gem_read_dist, k_min: params.gem_k_min,
         k_max: params.gem_k_max, k_seqs: params.gem_k_seqs,
         keep: keep, args: params.gem_args])

    // Peaks -> motifs
    TOP_PEAKS_FASTA(GEM.out.peaks, PREP_GENOME.out.fasta, PREP_GENOME.out.fai, PREP_GENOME.out.sizes,
                    [top_n: params.top_n_peaks, halfwidth: params.peak_halfwidth])
    def motif_db = params.memechip_db ? file(params.memechip_db, checkIfExists: true) : []
    MEME_CHIP(TOP_PEAKS_FASTA.out.fasta, motif_db, [args: params.memechip_args])
    MOTIF_QC(TOP_PEAKS_FASTA.out.windows.join(TOP_PEAKS_FASTA.out.ranked).join(MEME_CHIP.out.dir), motifs.files,
             [consensus: motifs.consensus, n_motifs: params.qc_motifs_n_denovo, pvalue: params.qc_motifs_pvalue,
              halfwidth: params.peak_halfwidth, top_n: params.top_n_peaks])

    // Peak -> gene annotation
    ANNOTATE_PEAKS(TOP_PEAKS_FASTA.out.ranked, PREP_GENOME.out.fasta, file(params.gtf, checkIfExists: true))

    // Coverage tracks (every library, controls too) and an IGV session over everything: each peak
    // set's ranked peaks, then the experiments' and the controls' bigWigs, in samplesheet order
    BIGWIG(SORT_FILTER_BAM.out.filtered,
           [extend: params.bigwig_extend, bin_size: params.bigwig_bin_size,
            normalize: params.bigwig_normalize, args: params.bigwig_args])
    def ctrl_ids   = rows.findAll { r -> r.control }*.control.toSet()
    // IGV can't read plain gzip (bgzip only with .fai + .gzi), so a .gz --fasta is replaced by
    // PREP_GENOME's decompressed copy, published to <outdir>/genome
    def fasta      = file(params.fasta)
    def gzipped    = fasta.name.endsWith('.gz')
    def igv_genome = gzipped ? file("${params.outdir}/genome/${fasta.baseName}").toAbsolutePath().toString()
                             : fasta.toRealPath().toString()
    if (!gzipped && !file("${igv_genome}.fai").exists())
        log.warn "IGV session: ${igv_genome} has no .fai next to it, so IGV can't load it (samtools faidx it)"
    IGV_SESSION([genome:  igv_genome,
                 gtf:     file(params.gtf).toRealPath().toString(),
                 peaks:   sets*.id,
                 bigwigs: rows.findAll { r -> !(r.sample in ctrl_ids) }*.sample + rows.findAll { r -> r.sample in ctrl_ids }*.sample])

    // QC report
    PEAK_QC(gem_in.map { meta, bams, bais, _cb, _ci -> [meta, bams, bais] }.join(TOP_PEAKS_FASTA.out.ranked))
    GEM_QC(GEM.out.dir)
    def qc = TRIM_GALORE.out.reports
        .mix(BOWTIE2_ALIGN.out.log, SORT_FILTER_BAM.out.qc, PEAK_QC.out, GEM_QC.out, MOTIF_QC.out.mqc,
             ANNOTATE_PEAKS.out.mqc)
        .flatten()
        .collect()
    MULTIQC(qc, file("${projectDir}/multiqc_config.yaml", checkIfExists: true))
}

// Samplesheet: sample,fastq_1[,fastq_2][,control][,group]
//   fastq_2 empty -> single-end. control = the `sample` id of another row (e.g. HaloTag/beads);
//   empty -> GEM runs without a control. Rows used as a control are aligned but not peak-called.
//   group = replicates of one experiment share a group name: besides each replicate's own peaks,
//   GEM is run once on all of them together (with all their controls), named after the group.
def read_samplesheet(path) {
    def raw = file(path, checkIfExists: true).splitCsv(header: true)
    if (!raw) error "Samplesheet is empty: ${path}"
    // A misspelt column (e.g. "Control") would otherwise be silently ignored
    def columns = ['sample', 'fastq_1', 'fastq_2', 'control', 'group']
    def unknown = raw[0].keySet() - columns
    if (unknown) error "Samplesheet: unknown column(s) ${unknown.join(', ')} (allowed: ${columns.join(', ')})"
    def missing = ['sample', 'fastq_1'] - raw[0].keySet()
    if (missing) error "Samplesheet: missing column(s) ${missing.join(', ')}"
    def rows = raw.collect { r ->
        [sample:  r.sample?.trim(),
         fastq_1: r.fastq_1?.trim(),
         fastq_2: r.fastq_2?.trim() ?: null,
         control: r.control?.trim() ?: null,
         group:   r.group?.trim() ?: null]
    }

    def ids = rows.collect { r -> r.sample }
    rows.each { r ->
        if (!r.sample || !(r.sample ==~ /[A-Za-z0-9._-]+/))
            error "Samplesheet: invalid sample id '${r.sample}' (allowed: letters, digits, . _ -)"
        if (!r.fastq_1)
            error "Samplesheet: sample '${r.sample}' has no fastq_1"
        if (ids.count(r.sample) > 1)
            error "Samplesheet: duplicate sample id '${r.sample}'"
        if (r.control && !(r.control in ids))
            error "Samplesheet: sample '${r.sample}' names control '${r.control}', which is not a sample"
        if (r.control == r.sample)
            error "Samplesheet: sample '${r.sample}' is its own control"
    }
    def ctrl_ids = rows.findAll { r -> r.control }.collect { r -> r.control }.toSet()
    rows.findAll { r -> r.sample in ctrl_ids && r.control }.each { r ->
        error "Samplesheet: '${r.sample}' is used as a control, so it can't have a control itself"
    }
    rows.findAll { r -> r.group }.groupBy { r -> r.group }.each { g, members ->
        // the group's peaks are published under its name, next to the samples'
        if (!(g ==~ /[A-Za-z0-9._-]+/))
            error "Samplesheet: invalid group '${g}' (allowed: letters, digits, . _ -)"
        if (g in ids)
            error "Samplesheet: group '${g}' has the same name as a sample"
        members.findAll { r -> r.sample in ctrl_ids }.each { r ->
            error "Samplesheet: '${r.sample}' is used as a control, so it can't be in a group"
        }
        if (members.any { r -> r.control } && !members.every { r -> r.control })
            error "Samplesheet: group '${g}' mixes samples with and without a control"
        if (members.size() == 1)
            log.warn "Samplesheet: group '${g}' has one sample, so its peaks repeat that sample's"
    }
    return rows
}

// Peak sets as [id, expt: [sample ids], ctrl: [sample ids]]: every sample that isn't a control,
// then every group (its samples and their distinct controls)
def peak_sets(rows) {
    def ctrl_ids = rows.findAll { r -> r.control }.collect { r -> r.control }.toSet()
    def expts    = rows.findAll { r -> !(r.sample in ctrl_ids) }
    def samples  = expts.collect { r -> [id: r.sample, expt: [r.sample], ctrl: r.control ? [r.control] : []] }
    def groups   = expts.findAll { r -> r.group }.groupBy { r -> r.group }.collect { g, members ->
        [id: g, expt: members*.sample, ctrl: members.findAll { r -> r.control }*.control.unique()]
    }
    return samples + groups
}

// --qc_motifs: IUPAC consensus strings and/or MEME-format motif files, as "CAGGTAG",
// "CAGGTAG,data/zelda.meme" or a YAML list -> [consensus: [strings], files: [paths]]
def qc_motifs_inputs() {
    def m     = params.qc_motifs
    def items = !m ? [] : m instanceof List ? m*.toString() : m.toString().tokenize(', ')
    def consensus = items.findAll { s -> s ==~ /(?i)[ACGTRYSWKMBDHVN]+/ }
    def files = (items - consensus).collect { s ->
        def f = file(s)
        if (!f.exists())
            error "--qc_motifs: '${s}' is neither a DNA consensus (IUPAC letters ACGT RYSWKM BDHV N) nor an existing motif file"
        f
    }
    if (files*.name.unique().size() < files.size())
        error "--qc_motifs: motif files need distinct file names"
    return [consensus: consensus*.toUpperCase(), files: files]
}

// --keep_chromosomes as a list; accepts "1,2,3", "1 2 3" or a YAML list. Empty list = all contigs.
def keep_chromosomes() {
    def k = params.keep_chromosomes
    if (!k) return []
    k instanceof List ? k*.toString() : k.toString().tokenize(', ')
}
