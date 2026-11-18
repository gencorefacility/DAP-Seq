// DAP-seq processes. Containers, resources and SLURM settings live in nextflow.config.
//
// Parameters reach each process as a `val opts` map built in main.nf, never as params.* in the
// script: Nextflow 26's parser leaves params out of the -resume cache key, so a changed param would
// silently reuse old results. Input values are always part of the key.

// ---------------------------------------------------------------------------
// Reference preparation
// ---------------------------------------------------------------------------

process PREP_GENOME {
    publishDir "${params.outdir}/genome", mode: 'copy', enabled: params.save_reference
    // A .gz input's decompressed copy + .fai are published regardless: the IGV session uses them,
    // since IGV can't read a gzipped FASTA (main.nf). The input is checked in saveAs, which runs per
    // task: enabled/pattern are evaluated once, without the inputs ("No such variable", lint passes).
    publishDir "${params.outdir}/genome", mode: 'copy', enabled: !params.save_reference,
        saveAs: { fn -> fasta.name.endsWith('.gz') && fn in [plain_name(fasta), plain_name(fasta) + '.fai'] ? fn : null }

    input:
    path fasta
    val opts      // keep: contigs to keep (empty = all)

    output:
    path "${plain_name(fasta)}", includeInputs: true, emit: fasta   // the input itself unless gzipped
    path "${plain_name(fasta)}.fai",                  emit: fai
    path "genome.sizes",                              emit: sizes
    path "chromosomes",                               emit: chr_dir

    script:
    // A gzipped FASTA (plain gzip or bgzip) is decompressed once here, since faidx and getfasta
    // can't read plain gzip; everything downstream takes this process's `fasta` output.
    // GEM's --genome wants a directory of per-sequence FASTAs, and GEM prefixes "chr" to contig
    // names when looking them up: contig "1" must be in chr1.fa. Names already starting with
    // "chr" are used as-is. GEM then reports "chr1"; TOP_PEAKS_FASTA maps it back to the genome's name.
    // --keep_chromosomes is checked here because samtools view silently skips unknown regions
    // (exit 0), so a typo would give empty BAMs. GEM strips a leading "chrom"/"chromosome" from
    // --g names but only "chr" from BAM names, so contigs named like that never match inside GEM.
    def keep  = opts.keep.join(' ')
    def plain = plain_name(fasta)
    def unzip = plain != fasta.name ? "gzip -cd ${fasta} > ${plain}" : ''
    """
    ${unzip}
    samtools faidx ${plain}
    cut -f1,2 ${plain}.fai > genome.sizes
    cut -f1 ${plain}.fai > contigs.txt

    # "1" and "chr1" would both be chr1.fa, and GEM couldn't tell them apart anyway
    dups=\$(sed 's/^chr//' contigs.txt | sort | uniq -d)
    [[ -z \$dups ]] || { echo "Contig names differ only by a leading 'chr', which GEM ignores:" \$dups >&2; exit 1; }

    mkdir chromosomes
    while read -r name; do
        [[ \$name == chr* ]] && file=\$name || file=chr\$name
        samtools faidx ${plain} "\$name" > "chromosomes/\$file.fa"
    done < contigs.txt

    if [[ -n "${keep}" ]]; then
        printf '%s\\n' ${keep} > used.txt
        missing=\$(grep -vxFf contigs.txt used.txt || true)
        [[ -z \$missing ]] || { echo "--keep_chromosomes not in ${fasta}:" \$missing >&2; exit 1; }
    else
        cp contigs.txt used.txt
    fi
    bad=\$(grep '^chrom' used.txt || true)
    [[ -z \$bad ]] || { echo "GEM can't handle contig names starting with 'chrom'; rename them in the FASTA or leave them out of --keep_chromosomes:" \$bad >&2; exit 1; }
    """
}

process BOWTIE2_BUILD {
    publishDir "${params.outdir}/genome", mode: 'copy', enabled: params.save_reference

    input:
    path fasta

    output:
    path "bowtie2_index"

    script:
    """
    mkdir bowtie2_index
    bowtie2-build --threads ${task.cpus} ${fasta} bowtie2_index/${fasta.baseName}
    """
}

// ---------------------------------------------------------------------------
// Reads -> filtered BAM
// ---------------------------------------------------------------------------

process TRIM_GALORE {
    tag "${meta.id}"
    publishDir "${params.outdir}/trim_galore", mode: 'copy', pattern: '*.{html,zip,txt}'

    input:
    tuple val(meta), path(reads)
    val opts      // args

    output:
    tuple val(meta), path("*.fq.gz"), emit: reads
    path "*.{html,zip,txt}",          emit: reports

    script:
    // Input goes to Trim Galore as plain in/<id>.fastq (<id>_1/_2 if paired): it names its reports
    // after the input file, and MultiQC takes the sample name from there. Gzipped input is
    // decompressed here (and removed after) because the pinned cutadapt 3.3 fails on some valid
    // BGZF files (Element AVITI) with a spurious "CRC check failed"; gzip reads them fine.
    // --gzip: gzipped output whatever the input.
    def paired = meta.single_end ? '' : '--paired'
    def names  = meta.single_end ? [meta.id] : ["${meta.id}_1", "${meta.id}_2"]
    def fqs    = [reads].flatten()
    def plain  = names.collect { n -> "in/${n}.fastq" }
    def prep   = fqs.indices.collect { i ->
        fqs[i].name.endsWith('.gz') ? "gzip -cd ${fqs[i]} > ${plain[i]}" : "ln -s ../${fqs[i]} ${plain[i]}"
    }.join('\n')
    """
    mkdir in
    ${prep}
    trim_galore -q 20 --fastqc --gzip --cores ${task.cpus} ${paired} --basename ${meta.id} ${opts.args} ${plain.join(' ')}
    rm -f ${plain.join(' ')}
    """
}

process BOWTIE2_ALIGN {
    tag "${meta.id}"
    publishDir "${params.outdir}/bowtie2", mode: 'copy', pattern: '*.log'

    input:
    tuple val(meta), path(reads)
    path index
    val opts      // args

    output:
    tuple val(meta), path("${meta.id}.sam"), emit: sam
    path "${meta.id}.bowtie2.log",           emit: log

    script:
    def input = meta.single_end ? "-U ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    INDEX=\$(find -L ${index}/ -name '*.rev.1.bt2*' | sed 's/\\.rev\\.1\\.bt2l\\?\$//')
    [[ \$(grep -c . <<< "\$INDEX") -eq 1 ]] || { echo "Expected one Bowtie2 index in ${index}, found: \${INDEX:-none}" >&2; exit 1; }
    bowtie2 -p ${task.cpus} -x \$INDEX ${input} ${opts.args} -S ${meta.id}.sam 2> ${meta.id}.bowtie2.log \\
        || { cat ${meta.id}.bowtie2.log >&2; exit 1; }
    """
}

process SORT_FILTER_BAM {
    tag "${meta.id}"
    publishDir "${params.outdir}/bam",      mode: 'copy', pattern: '*.{bam,bai}'
    publishDir "${params.outdir}/samtools", mode: 'copy', pattern: '*.{stats,flagstat,idxstats}'

    input:
    tuple val(meta), path(sam)
    val opts      // min_mapq, keep (empty = all contigs)

    output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"),                   emit: bam
    tuple val(meta), path("${meta.id}.filtered.bam"), path("${meta.id}.filtered.bam.bai"), emit: filtered
    path "*.{stats,flagstat,idxstats}",                                                    emit: qc

    script:
    // Upstream: MAPQ >= 30 on Arabidopsis nuclear chromosomes (1-5). Here both are params;
    // no keep_chromosomes = keep every contig.
    // -F 4 drops unmapped reads explicitly: -q only removes them because they have MAPQ 0, so with
    // --min_mapq 0 they'd reach GEM, FRiP and the bigWigs. A no-op at upstream's 30.
    // Paired-end: both mates of a fragment must pass (issue #1). Keep proper pairs (-f 2), drop unmapped,
    // mate-unmapped, secondary, QC-fail and duplicate records (-F 1804), MAPQ per record (-q). fixmate -r
    // then marks a record whose mate was dropped as mate-unmapped, and the second view removes it. The
    // filtered BAM keeps read 1 of each pair (-f 64), one record per DNA fragment: GEM has no paired-end
    // mode and would count both mates, i.e. every fragment twice. Read 1's record keeps the mate's
    // position and TLEN (the fragment for bigWigs).
    def regions = opts.keep.join(' ')
    """
    samtools view -bSu ${sam} | samtools sort -@ ${task.cpus} -o ${meta.id}.bam -
    samtools index ${meta.id}.bam
    if ${meta.single_end}; then
        samtools view -u -F 4 -q ${opts.min_mapq} ${meta.id}.bam ${regions} | samtools sort -@ ${task.cpus} -o ${meta.id}.filtered.bam -
    else
        samtools view -u -F 1804 -f 2 -q ${opts.min_mapq} ${meta.id}.bam ${regions} | samtools sort -n -@ ${task.cpus} -o pairs.bam -
        samtools fixmate -r pairs.bam pairs.fixmate.bam
        samtools view -u -F 1804 -f 2 pairs.fixmate.bam | samtools sort -@ ${task.cpus} -o pairs.filtered.bam -
        samtools view -b -f 64 -F 4 pairs.filtered.bam -o ${meta.id}.filtered.bam
    fi
    samtools index ${meta.id}.filtered.bam
    for b in ${meta.id} ${meta.id}.filtered; do
        samtools stats \$b.bam > \$b.stats
        samtools flagstat \$b.bam > \$b.flagstat
    done
    samtools idxstats ${meta.id}.bam > ${meta.id}.idxstats   # all contigs, incl. organelles
    """
}

// ---------------------------------------------------------------------------
// Peak calling. Every caller must emit a narrowPeak (col 10 = summit offset) so that
// everything downstream is caller-agnostic. Callers get SORT_FILTER_BAM's filtered BAMs, which for
// paired-end libraries hold read 1 only (with the mate's position and TLEN): enough for a
// fragment-based mode such as MACS3 -f BAMPE (it reads only read 1 anyway), not for a tool that
// needs both mates' records.
// ---------------------------------------------------------------------------

// One peak set = one sample, or a replicate group (samplesheet `group`): GEM's multi-replicate mode
// pools the replicates' reads (and their controls') and calls one set of events.
process GEM {
    tag "${meta.id}"
    publishDir "${params.outdir}/gem", mode: 'copy'

    input:
    tuple val(meta), path(bams), path(bais), path(ctrl_bams), path(ctrl_bais)   // lists; ctrl may be empty
    path sizes
    path chr_dir
    val opts      // jar, read_dist, k_min, k_max, k_seqs, keep, args

    output:
    tuple val(meta), path("${meta.id}/${meta.id}.GEM_events.narrowPeak"), emit: peaks
    tuple val(meta), path("${meta.id}"),                                   emit: dir

    script:
    // GEM's stdout is kept as <id>/<id>.GEM_log.txt: its read counts (before and after GEM's duplicate
    // filter) are only printed there, and GEM_QC reads them.
    // --expt<COND>/--ctrl<COND> must share the same condition suffix for GEM to pair them; repeating
    // either flag adds replicates to that condition (GEM pools repeated --ctrl too: verified 2026-10-22).
    // --t 1: GEM is single-threaded because more threads make its results non-reproducible (issue #2).
    def expt = [bams].flatten().collect { b -> "--expt${meta.id} ${b}" }.join(' ')
    def ctrl = [ctrl_bams].flatten().collect { b -> "--ctrl${meta.id} ${b}" }.join(' ')
    def subs = opts.keep ? "--subs \"${opts.keep.join(' ')}\"" : ''
    // Heap exhaustion would exit 1 (not retried); ExitOnOutOfMemoryError makes it exit 3, which
    // nextflow.config retries with more memory.
    def heap = (task.memory.toMega() * 0.85) as int   // leave headroom for JVM overhead
    """
    java -Xmx${heap}m -XX:+ExitOnOutOfMemoryError -jar ${opts.jar} \\
        --d ${opts.read_dist} \\
        --g ${sizes} \\
        --genome ${chr_dir} \\
        ${expt} \\
        ${ctrl} \\
        --f SAM \\
        --out ${meta.id} \\
        --t 1 \\
        --k_min ${opts.k_min} --k_max ${opts.k_max} --k_seqs ${opts.k_seqs} \\
        ${subs} \\
        --outNP --outBED --outMEME --outJASPAR --outHOMER \\
        --k_neg_dinu_shuffle --print_bound_seqs --print_aligned_seqs \\
        ${opts.args} | tee gem.log
    mv gem.log ${meta.id}/${meta.id}.GEM_log.txt
    """
}

// GEM's read counts before/after its duplicate filter, for MultiQC. Separate from GEM so changing the report never reruns peak calling.
process GEM_QC {
    tag "${meta.id}"

    input:
    tuple val(meta), path(gem)   // GEM's output dir

    output:
    path "${meta.id}.*_mqc.json"

    script:
    """
    gem_qc.py --sample ${meta.id} --dir ${gem}
    """
}

// ---------------------------------------------------------------------------
// Peaks -> motifs
// ---------------------------------------------------------------------------

process TOP_PEAKS_FASTA {
    tag "${meta.id}"
    publishDir "${params.outdir}/top_peaks", mode: 'copy'

    input:
    tuple val(meta), path(peaks)
    path fasta
    path fai
    path sizes
    val opts      // top_n (0 = all), halfwidth

    output:
    tuple val(meta), path("${meta.id}.ranked.narrowPeak"), emit: ranked    // all peaks, genome contig names
    tuple val(meta), path("${meta.id}.top.narrowPeak"),    emit: peaks
    tuple val(meta), path("${meta.id}.top.fa"),            emit: fasta
    tuple val(meta), path("${meta.id}.windows.fa"),        emit: windows   // all peaks, for MOTIF_QC

    script:
    // 1. Contig names: if a peak's contig isn't in the genome but is once a leading "chr" is
    //    dropped, use that; anything else is an error (never a blind strip).
    // 2. Rank by qValue (col 9) then signalValue (col 7), keep the top N (0 = all). Exact ties
    //    (common: GEM caps qValue at 999) are broken by the whole line as bytes (-k1, LC_ALL=C):
    //    that is GNU sort's implicit last resort, so the order matches upstream's
    //    `sort -k9,9nr -k7,7nr`, but explicit, so busybox sort and any host locale agree too.
    // 3. Re-centre each peak on its summit (start + col 10): [summit - hw, summit + hw + 1), so
    //    MEME-ChIP gets equal-width windows whatever the caller. GEM events are already exactly
    //    that with hw = 100 (201 bp, summit offset 100), so upstream's sequences are unchanged.
    // 4. slop -b 0 clamps coordinates to contig ends so getfasta doesn't fail on edge peaks.
    // 5. windows.fa: the same summit windows for every peak, for MOTIF_QC. Only full-width ones
    //    (peaks at contig ends are dropped) so match positions are comparable; header = rank.
    //    Done here so MOTIF_QC needs no bedtools, only MEME tools + Python.
    """
    awk -v FS='\\t' -v OFS='\\t' '
        NR == FNR { ok[\$1] = 1; next }
        !(\$1 in ok) {
            c = \$1; sub(/^chr/, "", c)
            if (c in ok) \$1 = c
            else { print "contig not in genome: " \$1 > "/dev/stderr"; exit 1 }
        }
        { print }' ${fai} ${peaks} \\
        | LC_ALL=C sort -k9,9nr -k7,7nr -k1 > ${meta.id}.ranked.narrowPeak

    awk -v FS='\\t' -v OFS='\\t' -v n=${opts.top_n} -v hw=${opts.halfwidth} '
        n && NR > n { exit }
        NF < 10 || \$10 < 0 { print "peak has no summit (col 10): " \$0 > "/dev/stderr"; exit 1 }
        { s = \$2 + \$10; \$2 = (s > hw ? s - hw : 0); \$3 = s + hw + 1; \$10 = s - \$2; print }' \\
        ${meta.id}.ranked.narrowPeak > ${meta.id}.top.narrowPeak

    bedtools slop -b 0 -i ${meta.id}.top.narrowPeak -g ${sizes} \\
        | bedtools getfasta -fi ${fasta} -bed stdin -fo ${meta.id}.top.fa

    awk -v FS='\\t' -v OFS='\\t' -v hw=${opts.halfwidth} '
        NR == FNR { len[\$1] = \$2; next }
        { s = \$2 + \$10 }
        \$10 >= 0 && s >= hw && s + hw + 1 <= len[\$1] { print \$1, s - hw, s + hw + 1, FNR }' \\
        ${sizes} ${meta.id}.ranked.narrowPeak \\
        | bedtools getfasta -fi ${fasta} -bed stdin -nameOnly -fo ${meta.id}.windows.fa
    """
}

// Peak count and FRiP (fraction of filtered reads overlapping a peak) for MultiQC. A replicate
// group's FRiP pools its replicates' reads, as GEM did.
process PEAK_QC {
    tag "${meta.id}"

    input:
    tuple val(meta), path(bams), path(bais), path(peaks)

    output:
    path "${meta.id}.peak_qc.tsv"

    script:
    """
    total=0; in_peaks=0
    for b in ${[bams].flatten().join(' ')}; do
        total=\$(( total + \$(samtools view -c \$b) ))
        if [[ -s ${peaks} ]]; then in_peaks=\$(( in_peaks + \$(samtools view -c -L ${peaks} \$b) )); fi
    done
    awk -v OFS='\\t' -v id=${meta.id} -v n=\$(wc -l < ${peaks}) -v a=\$in_peaks -v t=\$total 'BEGIN {
        print "Sample", "peaks", "frip"
        print id, n, (t ? sprintf("%.2f", 100 * a / t) : 0) }' > ${meta.id}.peak_qc.tsv
    """
}

process MEME_CHIP {
    tag "${meta.id}"
    publishDir "${params.outdir}/meme_chip", mode: 'copy'

    input:
    tuple val(meta), path(seqs)
    path motif_db
    val opts      // args

    output:
    tuple val(meta), path("${meta.id}"), emit: dir

    script:
    // Without a motif DB this is pure de novo discovery (no Tomtom comparison to known motifs).
    // TMPDIR: MEME is an MPI build, and OpenMPI fails if the node's TMPDIR isn't visible
    // inside the container. meme-chip exits 0 even when its steps fail, so
    // check its progress log (missing = failed too) and fail loudly instead of publishing a partial result.
    def db = motif_db ? "-db ${motif_db}" : ''
    """
    export TMPDIR=\$PWD/tmp && mkdir -p \$TMPDIR
    meme-chip -oc ${meta.id} ${db} -meme-p 1 ${opts.args} ${seqs}
    log=${meta.id}/progress_log.txt
    if [[ ! -s \$log ]] || grep -q 'status: [1-9]' \$log; then
        echo "meme-chip steps failed:" >&2
        grep 'status: [1-9]' \$log >&2 || echo "  no \$log" >&2
        exit 1
    fi
    """
}

// Motif QC over ALL peaks (MEME-ChIP only sees the top N): does the motif sit on the summit, and
// how far down the ranking do peaks still carry it? Motifs, in report order: --qc_motifs consensus
// strings (exact match, either strand), every motif in the --qc_motifs MEME files, then the top
// --qc_motifs_n_denovo of MEME-ChIP's combined.meme (its motifs sorted by significance). Matrices go
// through FIMO. Also MEME-ChIP's top de novo and known motif (its summary.tsv) as a MultiQC table.
process MOTIF_QC {
    tag "${meta.id}"
    publishDir "${params.outdir}/motif_qc", mode: 'copy'

    input:
    tuple val(meta), path(windows), path(ranked), path(memechip)
    path qc_files, stageAs: 'qc_motifs/*'   // --qc_motifs MEME files (may be none)
    val opts      // consensus (list), n_motifs, pvalue, halfwidth, top_n

    output:
    path "${meta.id}.*_mqc.json",                                          emit: mqc
    tuple val(meta), path("${meta.id}.motif_{summary,peaks}.tsv"),         emit: tables
    tuple val(meta), path("${meta.id}.fimo.tsv"),                          emit: fimo

    script:
    // windows = TOP_PEAKS_FASTA's summit windows for every peak (header = rank = line in ranked,
    // which gives the per-peak table its coordinates). The shuffled
    // windows (same dinucleotide content; fixed seed) give the false-positive rate.
    // FIMO --text streams every match (the default mode keeps only the best ~100,000). Consensus
    // strings are counted in motif_qc.py instead: FIMO's p-value cutoff doesn't suit them (an exact
    // 6-mer, or a 7-mer with one R/Y, can never reach p < 1e-4).
    // Each motif file is scanned on its own (fimo_<i>.*), so motif IDs may repeat across files.
    def consensus = opts.consensus.collect { c -> "--consensus ${c}" }.join(' ')
    def sources   = [qc_files].flatten().collect { f -> "${f} -1" } + ["${memechip}/combined.meme ${opts.n_motifs}"]
    def scans     = sources.indices.collect { i -> "scan ${sources[i]} fimo_${i + 1}" }.join('\n')
    def sets      = sources.indices.collect { i -> "--fimo fimo_${i + 1}" }.join(' ')
    """
    fasta-shuffle-letters -kmer 2 -seed 1 -dna ${windows} > shuffled.fa

    # scan <MEME file> <top n motifs, -1 = all> <output prefix>: FIMO on the real and shuffled windows
    scan() {
        awk -v OFS='\\t' -v n=\$2 '/^MOTIF/ && (n < 0 || ++k <= n) { print \$2, \$3 }' \$1 > \$3.motifs.tsv
        local sel=\$(awk '{ printf "--motif %s ", \$1 }' \$3.motifs.tsv)
        if [[ -z \$sel ]]; then   # --qc_motifs_n_denovo 0, or MEME-ChIP found nothing
            : > \$3.real.tsv; : > \$3.shuffled.tsv; return
        fi
        fimo --text --skip-matched-sequence --verbosity 1 --thresh ${opts.pvalue} \$sel \$1 ${windows} > \$3.real.tsv
        fimo --text --skip-matched-sequence --verbosity 1 --thresh ${opts.pvalue} \$sel \$1 shuffled.fa > \$3.shuffled.tsv
    }
    ${scans}
    awk 'FNR > 1 || NR == 1' fimo_*.real.tsv > ${meta.id}.fimo.tsv

    motif_qc.py --sample ${meta.id} ${consensus} ${sets} --memechip ${memechip} \\
        --windows ${windows} --shuffled-windows shuffled.fa --peaks ${ranked} \\
        --halfwidth ${opts.halfwidth} --top-n ${opts.top_n}
    """
}

// ---------------------------------------------------------------------------
// Peak -> gene annotation
// ---------------------------------------------------------------------------

// HOMER annotatePeaks on every peak's summit, with the run's own genome FASTA + GTF (no HOMER
// genome package, so any organism works). Output rows are in peak-rank order; PeakID = the rank
// (line number in <id>.ranked.narrowPeak). The annotation categories go to MultiQC.
process ANNOTATE_PEAKS {
    tag "${meta.id}"
    publishDir "${params.outdir}/annotation", mode: 'copy', pattern: '*.annotated.tsv'

    input:
    tuple val(meta), path(peaks)   // ranked narrowPeak
    path fasta
    path gtf

    output:
    tuple val(meta), path("${meta.id}.annotated.tsv"), emit: tsv
    path "${meta.id}.annotation_mqc.tsv",               emit: mqc

    script:
    // Contig names must match between peaks and GTF, else HOMER silently calls everything intergenic.
    def plain = plain_name(gtf)
    def unzip = plain != gtf.name ? "gzip -cd ${gtf} > ${plain}" : ''
    """
    ${unzip}
    awk -v OFS='\\t' '{ s = \$2 + \$10 + 1; print FNR, \$1, s, s, "+" }' ${peaks} > summits.txt
    if [[ -s summits.txt ]]; then
        shared=\$(awk -v FS='\\t' 'NR == FNR { g[\$1] = 1; next } \$2 in g { n++ } END { print n + 0 }' ${plain} summits.txt)
        [[ \$shared -gt 0 ]] || { echo "No peak is on a contig in ${gtf}: contig names differ from the genome's?" >&2; exit 1; }
        annotatePeaks.pl summits.txt ${fasta} -gtf ${plain} > annotated.txt
        # HOMER exits 0 even when its feature-assignment helper fails (every Annotation = NA)
        found=\$(awk -v FS='\\t' 'NR > 1 && \$8 != "NA" { n++ } END { print n + 0 }' annotated.txt)
        [[ \$found -gt 0 ]] || { echo "HOMER assigned no genomic features; see its log above" >&2; exit 1; }
        { head -n 1 annotated.txt; tail -n +2 annotated.txt | sort -t \$'\\t' -k1,1n; } > ${meta.id}.annotated.tsv
    else
        : > ${meta.id}.annotated.tsv
    fi

    # Peaks per category ("intron (FBtr..., intron 2 of 5)" -> intron) as a MultiQC bar graph
    awk -v FS='\\t' -v OFS='\\t' -v id=${meta.id} '
        FNR > 1 { c = \$8; sub(/ \\(.*/, "", c); n[c == "" ? "unannotated" : c]++ }
        END {
            print "# id: peak_annotation"
            print "# section_name: Peak annotation"
            print "# description: Genomic feature at each peak summit (HOMER annotatePeaks with the run GTF)."
            print "# plot_type: bargraph"
            print "# pconfig:"
            print "#     id: peak_annotation_plot"
            print "#     title: Peak annotation"
            print "#     ylab: Peaks"
            split("promoter-TSS exon intron TTS Intergenic", order, " ")   # what HOMER gives with a GTF
            hdr = "Sample"; row = id
            for (i = 1; i in order; i++) {
                c = order[i]; hdr = hdr OFS c; row = row OFS n[c] + 0; delete n[c] }
            other = 0; for (c in n) other += n[c]
            print hdr, "other"; print row, other
        }' ${meta.id}.annotated.tsv > ${meta.id}.annotation_mqc.tsv
    """
}

// ---------------------------------------------------------------------------
// Genome browser: bigWig tracks + IGV session
// ---------------------------------------------------------------------------

// Coverage from the filtered BAM (the reads GEM sees), with upstream's bamCoverage settings:
// --extendReads 200 --binSize 1 --normalizeUsing RPKM. Paired-end reads (read 1 only, see
// SORT_FILTER_BAM) are extended to their fragment (the mate's position, TLEN) instead of a fixed
// length. Upstream's --ignoreForNormalization Mt Pt isn't needed: the filtered BAM only has the
// --keep_chromosomes contigs.
process BIGWIG {
    tag "${meta.id}"
    publishDir "${params.outdir}/bigwig", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    val opts      // extend, bin_size, normalize, args

    output:
    path "${meta.id}.bw"

    script:
    // TMPDIR: the node's TMPDIR isn't visible inside the container
    def extend = meta.single_end ? "--extendReads ${opts.extend}" : '--extendReads'
    """
    export TMPDIR=\$PWD/tmp && mkdir -p \$TMPDIR
    bamCoverage --bam ${bam} --outFileName ${meta.id}.bw --outFileFormat bigwig \\
        --numberOfProcessors ${task.cpus} ${extend} --binSize ${opts.bin_size} \\
        --normalizeUsing ${opts.normalize} ${opts.args}
    """
}

// IGV desktop session (File > Open Session): track paths are relative to <outdir>/igv/, so the
// results directory can be copied elsewhere as a whole. The genome (--fasta) and genes (--gtf) are
// absolute paths: they only resolve where the reference lives. A .gz --fasta is replaced by its
// decompressed copy in <outdir>/genome, also by absolute path.
process IGV_SESSION {
    publishDir "${params.outdir}/igv", mode: 'copy'

    input:
    val opts      // genome, gtf (real paths), peaks (peak set ids), bigwigs (library ids)

    output:
    path "igv_session.xml"

    script:
    def tracks = opts.peaks.collect { id -> [path: "../top_peaks/${id}.ranked.narrowPeak", name: "${id} peaks"] } +
                 opts.bigwigs.collect { id -> [path: "../bigwig/${id}.bw", name: id] } +
                 [[path: opts.gtf, name: 'genes']]
    def resources = tracks.collect { t -> "        <Resource name=\"${t.name}\" path=\"${t.path}\"/>" }.join('\n')
    """
    cat <<'XML' > igv_session.xml
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="${opts.genome}" locus="All" version="8">
    <Resources>
${resources}
    </Resources>
</Session>
XML
    """
}

// ---------------------------------------------------------------------------
// Report
// ---------------------------------------------------------------------------

process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path files, stageAs: 'in/*'   // everything MultiQC should parse
    path config

    output:
    path "multiqc_report.html"
    path "multiqc_report_data"

    script:
    // -n: the config's title would otherwise rename the outputs (<title>_multiqc_report.html)
    """
    multiqc --config ${config} -n multiqc_report.html in/
    """
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

// genome.fa.gz -> genome.fa; anything else unchanged
def plain_name(file) {
    file.name.endsWith('.gz') ? file.baseName : file.name
}
