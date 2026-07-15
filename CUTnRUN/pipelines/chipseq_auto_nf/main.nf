#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def getRef(String species) {
    def ref = params.references[species]
    if (!ref) {
        error "[ERROR] Unsupported species: ${species}. Available keys: ${params.references.keySet().join(', ')}"
    }
    def required = ['bowtie2_index','gene_saf','te_saf','gene_anno','te_anno','blacklist']
    required.each { k ->
        if (!ref[k]) {
            error "[ERROR] Missing reference setting for species=${species}: ${k}"
        }
    }
    if (params.make_rpgc_track && !ref.effective_genome_size) {
        error "[ERROR] Missing reference setting for species=${species}: effective_genome_size"
    }
    return ref
}

workflow {
    log.info """
    ============================================================
            CHIP-seq / CUT&RUN core auto pipeline (Nextflow)
    ============================================================
    manifest   : ${params.manifest}
    outdir     : ${params.outdir}
    trim       : ${params.trim}
    analysis   : ${params.run_analysis}
    RPGC track : ${params.make_rpgc_track}
    TE align   : -k ${params.te_k}
    TE locus   : ${params.make_te_locus_best_track ? 'one-best-location CPM track' : 'disabled'}
    ============================================================
    """

    if (!params.manifest) {
        error "[ERROR] Please provide --manifest manifest.csv"
    }

    manifest_ch = Channel.fromPath(params.manifest, checkIfExists: true)

    RESOLVE_MANIFEST(manifest_ch)

    resolved_manifest_ch = RESOLVE_MANIFEST.out.resolved_manifest

    analysis_ref_ch = resolved_manifest_ch
        .splitCsv(header: true)
        .map { row -> row.species?.toString()?.trim() }
        .first()
        .map { species -> getRef(species) }

    sample_reads_ch = resolved_manifest_ch
        .splitCsv(header: true)
        .map { row ->
            def sample = row.sample?.toString()?.trim()
            if (!sample) {
                error "[ERROR] sample column is required in manifest"
            }
            def species = row.species?.toString()?.trim()
            if (!species) {
                error "[ERROR] species column is required for sample=${sample}"
            }
            def ref = getRef(species)
            def fastq1 = row.fastq_1?.toString()?.trim()
            def fastq2 = row.fastq_2?.toString()?.trim()
            if (!fastq1) {
                error "[ERROR] fastq_1 is empty after resolving manifest for sample=${sample}"
            }
            def reads = fastq2 ? [ file(fastq1, checkIfExists: true), file(fastq2, checkIfExists: true) ]
                              : [ file(fastq1, checkIfExists: true) ]
            def layout = row.layout?.toString()?.trim()?.toUpperCase()
            if (!layout) {
                layout = fastq2 ? 'PE' : 'SE'
            }
            def isIgg = (row.is_igg?.toString()?.trim()?.toLowerCase() in ['1','true','yes','y'])
            def meta = [
                sample    : sample,
                species   : species,
                assay     : (row.assay?.toString()?.trim() ?: params.assay),
                group     : (row.group?.toString()?.trim() ?: sample),
                replicate : (row.replicate?.toString()?.trim() ?: ''),
                igg       : (row.igg?.toString()?.trim() ?: ''),
                is_igg    : isIgg,
                layout    : layout,
                paired    : layout == 'PE',
                ref       : ref
            ]
            tuple(meta, reads)
        }

    PREP_READS(sample_reads_ch)

    prepped_reads_ch = PREP_READS.out.prepped_reads

    ALIGN(prepped_reads_ch)
    SORT_BAM(ALIGN.out.aligned_bam)
    CLEAN_BAM(SORT_BAM.out.sorted_bam)
    clean_bam_ch = CLEAN_BAM.out.clean_bam

    ALIGN_TE(prepped_reads_ch)
    SORT_BAM_TE(ALIGN_TE.out.aligned_bam)
    CLEAN_BAM_TE(SORT_BAM_TE.out.sorted_bam)
    te_bam_ch = CLEAN_BAM_TE.out.te_bam

    if (params.make_te_locus_best_track) {
        TE_LOCUS_BEST_BAM(SORT_BAM_TE.out.sorted_bam)
        BIGWIG_TE_LOCUS_BEST(TE_LOCUS_BEST_BAM.out.locus_bam)
    }

    BIGWIG(clean_bam_ch)

    if (params.make_te_tracks) {
        BIGWIG_TE(te_bam_ch)
    }

    peak_target_ch = clean_bam_ch
        .filter { meta, bam, bai -> !meta.is_igg }

    CALL_PEAKS(peak_target_ch)

    FEATURECOUNTS_GENE(clean_bam_ch.collect(flat: false))
    FEATURECOUNTS_TE(te_bam_ch.collect(flat: false))

    if (params.run_analysis) {
        RUN_ANALYSIS(
            FEATURECOUNTS_GENE.out.counts,
            FEATURECOUNTS_TE.out.counts,
            resolved_manifest_ch,
            analysis_ref_ch
        )
    }
}

process RESOLVE_MANIFEST {
    tag 'resolve_manifest'
    publishDir "${params.outdir}/manifest", mode: 'copy'

    input:
    path manifest

    output:
    path 'resolved_manifest.csv', emit: resolved_manifest

    script:
    """
    set -euo pipefail
    echo "[RESOLVE_MANIFEST] resolving manifest: ${manifest}"
    python ${projectDir}/bin/resolve_manifest.py \
      --manifest ${manifest} \
      --outdir ${params.outdir}/downloads \
      --threads ${params.download_threads} \
      --output resolved_manifest.csv
    echo "[RESOLVE_MANIFEST] done"
    """
}

process PREP_READS {
    tag { meta.sample }
    publishDir "${params.outdir}/01_prepared_reads", mode: 'copy', pattern: '*.fastq.gz'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.fastq.gz'), emit: prepped_reads

    script:
    def sample = meta.sample
    def layout = meta.layout
    def trimFlag = params.trim ? 'true' : 'false'
    def read1 = reads[0]
    def read2 = reads.size() > 1 ? reads[1] : null
    if (layout == 'PE' && !read2) {
        error "[ERROR] Sample ${sample} is marked PE but fastq_2 is missing"
    }
    if (layout == 'PE') {
        return """
        set -euo pipefail
        echo "[PREP_READS] sample=${sample} layout=${layout} trim=${trimFlag}"
        if [[ "${trimFlag}" == "true" ]]; then
            trim_galore \
              --paired \
              --gzip \
              --cores ${task.cpus} \
              --fastqc \
              ${read1} ${read2} \
              -o . \
              2> ${sample}.trim.log
            mv *_val_1*.fq.gz ${sample}_R1.trim.fastq.gz
            mv *_val_2*.fq.gz ${sample}_R2.trim.fastq.gz
        else
            ln -s ${read1} ${sample}_R1.trim.fastq.gz
            ln -s ${read2} ${sample}_R2.trim.fastq.gz
            echo "[PREP_READS] trim skipped for ${sample}" > ${sample}.trim.log
        fi
        """
    }
    else {
        return """
        set -euo pipefail
        echo "[PREP_READS] sample=${sample} layout=${layout} trim=${trimFlag}"
        if [[ "${trimFlag}" == "true" ]]; then
            trim_galore \
              --gzip \
              --cores ${task.cpus} \
              --fastqc \
              ${read1} \
              -o . \
              2> ${sample}.trim.log
            if ls *_trimmed.f*q.gz >/dev/null 2>&1; then
              mv *_trimmed.f*q.gz ${sample}.trim.fastq.gz
            else
              mv *_trimmed.fastq.gz ${sample}.trim.fastq.gz
            fi
        else
            ln -s ${read1} ${sample}.trim.fastq.gz
            echo "[PREP_READS] trim skipped for ${sample}" > ${sample}.trim.log
        fi
        """
    }
}

process ALIGN {
    tag { meta.sample }
    publishDir "${params.outdir}/02_align", mode: 'copy', pattern: '*.{bam,log}'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.sample}.aligned.bam"), emit: aligned_bam

    script:
    def sample = meta.sample
    def idx = meta.ref.bowtie2_index
    if (meta.paired) {
        return """
        set -euo pipefail
        echo "[ALIGN] sample=${sample} species=${meta.species} layout=PE index=${idx}"
        bowtie2 \
          -p ${task.cpus} \
          -t \
          -q \
          -N 1 \
          -L 25 \
          -X ${params.max_frag} \
          --rg-id ${sample} \
          --rg SM:${sample} \
          --no-unal \
          -x ${idx} \
          -1 ${reads[0]} \
          -2 ${reads[1]} \
          2> ${sample}.align.log \
        | samtools view -bS -o ${sample}.aligned.bam -
        echo "[ALIGN] done ${sample}"
        """
    } else {
        return """
        set -euo pipefail
        echo "[ALIGN] sample=${sample} species=${meta.species} layout=SE index=${idx}"
        bowtie2 \
          -p ${task.cpus} \
          -t \
          -q \
          -N 1 \
          -L 25 \
          --rg-id ${sample} \
          --rg SM:${sample} \
          --no-unal \
          -x ${idx} \
          -U ${reads[0]} \
          2> ${sample}.align.log \
        | samtools view -bS -o ${sample}.aligned.bam -
        echo "[ALIGN] done ${sample}"
        """
    }
}

process ALIGN_TE {
    tag { meta.sample }
    publishDir "${params.outdir}/02_align_te", mode: 'copy', pattern: '*.{bam,log}'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.sample}.aligned_te.bam"), emit: aligned_bam

    script:
    def sample = meta.sample
    def idx = meta.ref.bowtie2_index
    if (meta.paired) {
        return """
        set -euo pipefail
        echo "[ALIGN_TE] sample=${sample} layout=PE k=${params.te_k}"
        bowtie2 \
          -p ${task.cpus} \
          --very-sensitive \
          -t -q \
          -k ${params.te_k} \
          -X ${params.max_frag} \
          --rg-id ${sample} \
          --rg SM:${sample} \
          --no-unal \
          --no-mixed \
          --no-discordant \
          -x ${idx} \
          -1 ${reads[0]} \
          -2 ${reads[1]} \
          2> ${sample}.align_te.log \
        | samtools view -bS -o ${sample}.aligned_te.bam -
        """
    }
    return """
    set -euo pipefail
    echo "[ALIGN_TE] sample=${sample} layout=SE k=${params.te_k}"
    bowtie2 \
      -p ${task.cpus} \
      --very-sensitive \
      -t -q \
      -k ${params.te_k} \
      --rg-id ${sample} \
      --rg SM:${sample} \
      --no-unal \
      -x ${idx} \
      -U ${reads[0]} \
      2> ${sample}.align_te.log \
    | samtools view -bS -o ${sample}.aligned_te.bam -
    """
}

process SORT_BAM {
    tag { meta.sample }
    publishDir "${params.outdir}/03_sorted_bam", mode: 'copy', pattern: '*.{bam,bai}'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.sample}.sorted.bam"), path("${meta.sample}.sorted.bam.bai"), emit: sorted_bam

    script:
    def sample = meta.sample
    """
    set -euo pipefail
    echo "[SORT_BAM] sorting ${sample}"
    samtools sort -@ ${task.cpus} -o ${sample}.sorted.bam ${bam}
    samtools index ${sample}.sorted.bam
    echo "[SORT_BAM] done ${sample}"
    """
}

process SORT_BAM_TE {
    tag { meta.sample }
    publishDir "${params.outdir}/03_sorted_bam_te", mode: 'copy', pattern: '*.{bam,bai}'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.sample}.sorted_te.bam"), path("${meta.sample}.sorted_te.bam.bai"), emit: sorted_bam

    script:
    """
    set -euo pipefail
    samtools sort -@ ${task.cpus} -o ${meta.sample}.sorted_te.bam ${bam}
    samtools index ${meta.sample}.sorted_te.bam
    """
}

process CLEAN_BAM {
    tag { meta.sample }
    publishDir "${params.outdir}/04_clean_bam", mode: 'copy', pattern: '*'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.sample}_clean.bam"), path("${meta.sample}_clean.bam.bai"), emit: clean_bam

    script:
    def sample = meta.sample
    def blacklist = meta.ref.blacklist ?: ''
    def cleanRegex = meta.ref.clean_regex ?: ''
    def mapqCmd = meta.paired ? "samtools view -@ ${task.cpus} -h -b -q ${params.min_mapq} -F 1804 -f 2" : "samtools view -@ ${task.cpus} -h -b -q ${params.min_mapq} -F 1804"

    def cleanStep = cleanRegex ? """
    samtools view -h ${sample}.filtered.bam \\
      | awk '\$0 ~ /^@/ || \$3 !~ /${cleanRegex}/' \\
      | grep -v 'XS:i:' \\
      | samtools view -b - \\
      | samtools sort -@ ${task.cpus} -o ${sample}_clean.bam -
    """ : """
    samtools view -h ${sample}.filtered.bam \\
      | grep -v 'XS:i:' \\
      | samtools view -b - \\
      | samtools sort -@ ${task.cpus} -o ${sample}_clean.bam -
    """

    def blacklistStep = blacklist ? """
    echo "[CLEAN_BAM] removing blacklist for ${sample}: ${blacklist}"
    bedtools intersect -v -abam ${sample}.rmdup.bam -b ${blacklist} > ${sample}.blrm.bam
    input_bam=${sample}.blrm.bam
    """ : """
    echo "[CLEAN_BAM] blacklist skipped for ${sample}"
    cp ${sample}.rmdup.bam ${sample}.blrm.bam
    input_bam=${sample}.blrm.bam
    """

    """
    set -euo pipefail
    echo "[CLEAN_BAM] sample=${sample} layout=${meta.layout}"

    if samtools view -H ${bam} | grep -q '^@RG'; then
      markdup_input=${bam}
    else
      echo "[CLEAN_BAM] no @RG found; adding read group for ${sample}"
      samtools addreplacerg -r "ID:${sample}" -r "SM:${sample}" -o ${sample}.rg.bam ${bam}
      samtools index ${sample}.rg.bam
      markdup_input=${sample}.rg.bam
    fi

    picard MarkDuplicates \
      I=\$markdup_input \
      O=${sample}.rmdup.bam \
      M=${sample}.dup_metrics.txt \
      REMOVE_DUPLICATES=true \
      ASSUME_SORTED=true \
      CREATE_INDEX=true \
      VALIDATION_STRINGENCY=LENIENT \
      > ${sample}.markdup.log 2>&1

    ${blacklistStep}

    echo "[CLEAN_BAM] MAPQ/proper-pair filtering for ${sample}"
    ${mapqCmd} \$input_bam > ${sample}.filtered.bam

    echo "[CLEAN_BAM] removing unwanted contigs / XS tags for ${sample}"
    ${cleanStep}

    samtools index ${sample}_clean.bam
    samtools flagstat ${sample}_clean.bam > ${sample}_clean.flagstat.txt
    echo "[CLEAN_BAM] done ${sample}"
    """
}

process CLEAN_BAM_TE {
    tag { meta.sample }
    publishDir "${params.outdir}/04_te_bam", mode: 'copy', pattern: '*'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.sample}_te.bam"), path("${meta.sample}_te.bam.bai"), emit: te_bam

    script:
    def sample = meta.sample
    def blacklist = meta.ref.blacklist ?: ''
    def cleanRegex = meta.ref.clean_regex ?: ''
    def teMinMapq = params.te_min_mapq ?: 0
    def teExcludeFlags = params.te_exclude_flags ?: 4
    def teProperPair = params.te_proper_pair_only == null ? true : params.te_proper_pair_only.toString().toBoolean()
    def teRemoveBlacklist = params.te_remove_blacklist == null ? false : params.te_remove_blacklist.toString().toBoolean()
    def duplicatePolicy = (params.te_duplicate_policy ?: 'mark').toString().toLowerCase()
    if (!(duplicatePolicy in ['mark', 'keep', 'remove'])) {
        error "[CLEAN_BAM_TE] te_duplicate_policy must be mark, keep, or remove"
    }
    def pairArg = (teProperPair && meta.paired) ? "-f 2" : ""
    def blacklistStep = (teRemoveBlacklist && blacklist) ? """
    echo "[CLEAN_BAM_TE] removing blacklist for ${sample}: ${blacklist}"
    bedtools intersect -v -abam ${bam} -b ${blacklist} > ${sample}.te.blrm.bam
    """ : """
    echo "[CLEAN_BAM_TE] blacklist skipped for ${sample}"
    cp ${bam} ${sample}.te.blrm.bam
    """
    def duplicateStep
    if (duplicatePolicy == 'keep') {
        duplicateStep = """
        echo "[CLEAN_BAM_TE] duplicate policy=keep"
        cp ${sample}.te.queryname.bam ${sample}.te.dup_processed.bam
        """
    } else {
        def removeDuplicates = duplicatePolicy == 'remove' ? 'true' : 'false'
        duplicateStep = """
        echo "[CLEAN_BAM_TE] duplicate policy=${duplicatePolicy}"
        picard MarkDuplicates \
          I=${sample}.te.queryname.bam \
          O=${sample}.te.dup_processed.bam \
          M=${sample}.te.dup_metrics.txt \
          REMOVE_DUPLICATES=${removeDuplicates} \
          ASSUME_SORT_ORDER=queryname \
          CREATE_INDEX=false \
          VALIDATION_STRINGENCY=LENIENT \
          > ${sample}.te.markdup.log 2>&1
        """
    }

    """
    set -euo pipefail
    echo "[CLEAN_BAM_TE] sample=${sample} layout=${meta.layout}"
    echo "[CLEAN_BAM_TE] min_mapq=${teMinMapq} exclude_flags=${teExcludeFlags} proper_pair_only=${teProperPair}"

    ${blacklistStep}

    echo "[CLEAN_BAM_TE] relaxed filtering while retaining secondary alignments"
    samtools view \\
      -@ ${task.cpus} \\
      -h -b \\
      -q ${teMinMapq} \\
      -F ${teExcludeFlags} \\
      ${pairArg} \\
      ${sample}.te.blrm.bam \\
      > ${sample}.te.filtered.bam

    if [[ -n "${cleanRegex}" ]]; then
      samtools view -h ${sample}.te.filtered.bam \\
        | awk '\$0 ~ /^@/ || \$3 !~ /${cleanRegex}/' \\
        | samtools view -b - \\
        > ${sample}.te.contig.bam
    else
      cp ${sample}.te.filtered.bam ${sample}.te.contig.bam
    fi

    samtools sort -n -@ ${task.cpus} -o ${sample}.te.queryname.bam ${sample}.te.contig.bam
    ${duplicateStep}

    samtools view -h ${sample}.te.dup_processed.bam \\
      | python ${projectDir}/bin/add_nh_tag.py \\
      | samtools view -b - \\
      | samtools sort -@ ${task.cpus} -o ${sample}_te.bam -

    samtools index ${sample}_te.bam
    samtools flagstat ${sample}_te.bam > ${sample}_te.flagstat.txt
    nh_tagged="\$(samtools view -c -d NH ${sample}_te.bam)"
    if [[ "\${nh_tagged}" -eq 0 ]]; then
      echo "[CLEAN_BAM_TE] ERROR: no NH tags found in ${sample}_te.bam" >&2
      exit 1
    fi
    {
      printf 'metric\\tcount\\n'
      printf 'input\\t%s\\n' "\$(samtools view -c ${bam})"
      printf 'input_secondary\\t%s\\n' "\$(samtools view -c -f 256 ${bam})"
      printf 'final\\t%s\\n' "\$(samtools view -c ${sample}_te.bam)"
      printf 'final_secondary\\t%s\\n' "\$(samtools view -c -f 256 ${sample}_te.bam)"
      printf 'final_duplicate\\t%s\\n' "\$(samtools view -c -f 1024 ${sample}_te.bam)"
      printf 'final_nh_tagged\\t%s\\n' "\${nh_tagged}"
    } > ${sample}_te.clean_counts.tsv
    echo "[CLEAN_BAM_TE] done ${sample}"
    """
}

process BIGWIG {
    tag { meta.sample }
    publishDir "${params.outdir}/05_tracks", mode: 'copy', pattern: '*.bw'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    path '*.bw'

    script:
    def sample = meta.sample
    def extend = meta.paired ? '' : "--extendReads ${params.extend_reads_se}"
    def rpgcCmd = params.make_rpgc_track ? """
    bamCoverage \\
      --binSize 10 \\
      -p ${task.cpus} \\
      --normalizeUsing RPGC \\
      --effectiveGenomeSize ${meta.ref.effective_genome_size} \\
      --ignoreForNormalization ${params.ignore_for_normalization} \\
      -b ${bam} \\
      -o ${sample}_10bp_rpgc.bw
    """ : 'echo "[BIGWIG] RPGC track skipped"'

    """
    set -euo pipefail
    echo "[BIGWIG] sample=${sample} layout=${meta.layout}"
    bamCoverage \
      -b ${bam} \
      -o ${sample}_100bp_rpkm.bw \
      -bs 100 \
      -p ${task.cpus} \
      --normalizeUsing RPKM \
      ${extend}
    ${rpgcCmd}
    echo "[BIGWIG] done ${sample}"
    """
}

process BIGWIG_TE {
    tag { meta.sample }
    publishDir "${params.outdir}/05_tracks_te", mode: 'copy', pattern: '*.bw'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    path '*.bw'

    script:
    def sample = meta.sample
    def binSize = params.te_bw_binsize ?: 10
    def norm = (params.te_track_normalization ?: 'CPM').toString().toUpperCase()
    def extend = meta.paired ? '' : "--extendReads ${params.extend_reads_se}"

    def normArgs
    if (norm == 'RPGC') {
        if (!meta.ref.effective_genome_size) {
            error "[BIGWIG_TE] effective_genome_size is required when te_track_normalization=RPGC"
        }
        normArgs = "--normalizeUsing RPGC --effectiveGenomeSize ${meta.ref.effective_genome_size} --ignoreForNormalization ${params.ignore_for_normalization}"
    } else {
        normArgs = "--normalizeUsing ${norm}"
    }

    """
    set -euo pipefail
    echo "[BIGWIG_TE] sample=${sample} layout=${meta.layout} binSize=${binSize} normalization=${norm}"
    bamCoverage \\
      -b ${bam} \\
      -o ${sample}_te_${binSize}bp_${norm.toLowerCase()}.bw \\
      --binSize ${binSize} \\
      -p ${task.cpus} \\
      ${normArgs} \\
      ${extend}
    echo "[BIGWIG_TE] done ${sample}"
    """
}

process TE_LOCUS_BEST_BAM {
    tag { meta.sample }
    publishDir "${params.outdir}/04_te_locus_best_bam", mode: 'copy', pattern: '*'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.sample}_te_locus_best.bam"), path("${meta.sample}_te_locus_best.bam.bai"), emit: locus_bam

    script:
    def sample = meta.sample
    def properPair = meta.paired ? '-f 2' : ''
    """
    set -euo pipefail
    echo "[TE_LOCUS_BEST_BAM] sample=${sample} retaining one primary best alignment"
    samtools view \
      -@ ${task.cpus} \
      -h -b \
      -F 2308 \
      ${properPair} \
      ${bam} \
      > ${sample}_te_locus_best.bam
    samtools index ${sample}_te_locus_best.bam
    {
      printf 'metric\\tcount\\n'
      printf 'primary_best\\t%s\\n' "\$(samtools view -c ${sample}_te_locus_best.bam)"
      printf 'xs_tagged\\t%s\\n' "\$(samtools view ${sample}_te_locus_best.bam | grep -c "XS:i:" || true)"
    } > ${sample}_te_locus_best.metrics.tsv
    """
}

process BIGWIG_TE_LOCUS_BEST {
    tag { meta.sample }
    publishDir "${params.outdir}/05_tracks_te_locus_best", mode: 'copy', pattern: '*.bw'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    path '*.bw'

    script:
    def sample = meta.sample
    def binSize = params.te_locus_best_binsize ?: 5
    def fragmentArgs = meta.paired
        ? '--extendReads --samFlagInclude 64'
        : "--extendReads ${params.extend_reads_se}"
    """
    set -euo pipefail
    echo "[BIGWIG_TE_LOCUS_BEST] sample=${sample} binSize=${binSize} normalization=CPM"
    bamCoverage \
      -b ${bam} \
      -o ${sample}_te_locus_best_${binSize}bp_cpm.bw \
      --binSize ${binSize} \
      --normalizeUsing CPM \
      --exactScaling \
      -p ${task.cpus} \
      ${fragmentArgs}
    """
}

process CALL_PEAKS {
    tag { meta.sample }
    publishDir { "${params.outdir}/06_peaks/${meta.sample}" }, mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    path '*'

    script:
    def sample = meta.sample
    def peakFormat = meta.paired ? 'BAMPE' : 'BAM'
    def macsGenome = meta.ref.macs2_genome ?: (meta.species == 'hg38' ? 'hs' : 'mm')
    def peakCaller = params.macs_cmd ?: 'macs3'
    def macsCutoff = params.macs_pvalue ? "-p ${params.macs_pvalue}" : (params.macs_qvalue ? "-q ${params.macs_qvalue}" : "")
    def broadCutoff = params.macs_broad_cutoff ? "--broad-cutoff ${params.macs_broad_cutoff}" : ""
    def cutoffAnalysis = params.macs_cutoff_analysis ? '--cutoff-analysis' : ''
    // CLEAN_BAM publishes deterministic sample-named BAMs.  Referencing the
    // published control directly avoids a queue/value-channel pairing trap
    // when a single control map is reused for several target samples.
    def controlDir = "${params.outdir}/04_clean_bam"
    """
    set -euo pipefail
    echo "[CALL_PEAKS] sample=${sample} control=${meta.igg ?: 'NONE'} format=${peakFormat}"
    control_arg=""
    if [[ -n "${meta.igg}" ]]; then
      control_bam="${controlDir}/${meta.igg}_clean.bam"
      if [[ -s "\${control_bam}" ]]; then
        control_arg="-c \${control_bam}"
        echo "[CALL_PEAKS] using control bam: \${control_bam}"
      else
        echo "[CALL_PEAKS] warning: control BAM not found at \${control_bam}, continue without control"
      fi
    fi

    mkdir -p narrow broad

    ${peakCaller} callpeak \
      -t ${bam} \
      \${control_arg} \
      -f ${peakFormat} \
      -n ${sample} \
      -g ${macsGenome} \
      --keep-dup 1 \
      ${cutoffAnalysis} \
      ${macsCutoff} \
      --outdir narrow \
      > ${sample}.macs2.narrow.log 2>&1

    ${peakCaller} callpeak \
      -t ${bam} \
      \${control_arg} \
      -f ${peakFormat} \
      -n ${sample} \
      -g ${macsGenome} \
      --keep-dup 1 \
      ${cutoffAnalysis} \
      ${macsCutoff} \
      --broad \
      ${broadCutoff} \
      --outdir broad \
      > ${sample}.macs2.broad.log 2>&1

    [[ -f narrow/${sample}_peaks.narrowPeak ]] && cut -f1-3 narrow/${sample}_peaks.narrowPeak > narrow/${sample}_peaks.narrowPeak.bed || true
    [[ -f broad/${sample}_peaks.broadPeak ]] && cut -f1-3 broad/${sample}_peaks.broadPeak > broad/${sample}_peaks.broadPeak.bed || true
    [[ -f broad/${sample}_peaks.gappedPeak ]] && cut -f1-3 broad/${sample}_peaks.gappedPeak > broad/${sample}_peaks.gappedPeak.bed || true
    echo "[CALL_PEAKS] done ${sample}"
    """
}

process FEATURECOUNTS_GENE {
    tag 'featurecounts_gene'
    publishDir "${params.outdir}/07_featurecounts", mode: 'copy'

    input:
    val entries

    output:
    path 'featurecounts_gene.txt', emit: counts
    path 'featurecounts_gene.txt.summary'

    script:
    def firstMeta = entries[0][0]
    def saf = firstMeta.ref.gene_saf
    def bamList = entries.collect { it[1].toString() }.join(' ')
    def anyPaired = entries.any { it[0].paired }
    def allPaired = entries.every { it[0].paired }
    if (anyPaired && !allPaired) {
        error "[ERROR] Mixed PE/SE samples are not supported"
    }
    def pairedArg = allPaired ? '-p --countReadPairs' : ''
    """
    set -euo pipefail
    echo "[FEATURECOUNTS_GENE] annotation=${saf}"
    featureCounts \
      -a ${saf} \
      -o featurecounts_gene.txt \
      -F SAF \
      --ignoreDup \
      ${pairedArg} \
      -T ${task.cpus} \
      ${bamList}
    """
}

process FEATURECOUNTS_TE {
    tag 'featurecounts_te'
    publishDir "${params.outdir}/07_featurecounts", mode: 'copy'

    input:
    val entries

    output:
    path 'featurecounts_te.txt', emit: counts
    path 'featurecounts_te.txt.summary'

    script:
    def firstMeta = entries[0][0]
    def saf = firstMeta.ref.te_saf
    def bamList = entries.collect { it[1].toString() }.join(' ')
    def anyPaired = entries.any { it[0].paired }
    def allPaired = entries.every { it[0].paired }
    if (anyPaired && !allPaired) {
        error "[ERROR] Mixed PE/SE samples are not supported"
    }
    def pairedArg = allPaired ? '-p --countReadPairs' : ''
    """
    set -euo pipefail
    echo "[FEATURECOUNTS_TE] annotation=${saf}"
    featureCounts \
      -a ${saf} \
      -o featurecounts_te.txt \
      -F SAF \
      --ignoreDup \
      ${pairedArg} \
      -T ${task.cpus} \
      -M \
      --fraction \
      ${bamList}
    if awk 'BEGIN{FS="\\t"} !/^#/ && NR>2 {for(i=7;i<=NF;i++) if(\$i != int(\$i)){found=1; exit}} END{exit !found}' featurecounts_te.txt; then
      echo "[FEATURECOUNTS_TE] fractional counts detected"
    else
      echo "[FEATURECOUNTS_TE] warning: no fractional values detected" >&2
    fi
    """
}

process RUN_ANALYSIS {
    tag 'unified_analysis'
    publishDir "${params.outdir}/08_analysis", mode: 'copy'

    input:
    path gene_counts_file
    path te_counts_file
    path manifest_file
    val ref

    output:
    path 'results'

    script:
    """
    set -euo pipefail
    echo "[RUN_ANALYSIS] preparing config"
    python ${projectDir}/bin/build_analysis_config.py \
      --manifest ${manifest_file} \
      --gene-counts ${gene_counts_file} \
      --te-counts ${te_counts_file} \
      --output config.auto.R \
      --results-dir results \
      --repeat-anno ${ref.te_anno} \
      --gene-anno ${ref.gene_anno} \
      --deny-list ${ref.blacklist} \
      --gene-saf ${ref.gene_saf} \
      --te-saf ${ref.te_saf} \
      --te-classes ${params.te_classes.join(',')} \
      --min-total-normalized-counts ${params.min_total_normalized_counts} \
      --n-labels-genes ${params.n_labels_genes} \
      --n-labels-te-family ${params.n_labels_te_family} \
      --n-labels-te-repname ${params.n_labels_te_repname} \
      --max-te-plot-rows ${params.max_te_plot_rows} \
      --max-te-scatter-points ${params.max_te_scatter_points} \
      --max-te-facet-levels ${params.max_te_facet_levels} \
      --analysis-code ${projectDir}/analysis/run.R \
      --analysis-code ${projectDir}/analysis/scripts/01_prepare_data.R \
      --analysis-code ${projectDir}/analysis/scripts/02_normalize_and_process.R \
      --analysis-code ${projectDir}/../count_draw/scripts/03_generate_visuals.R \
      --fingerprint-param te_classes=${params.te_classes.join(',')} \
      --fingerprint-param min_total_normalized_counts=${params.min_total_normalized_counts}

    Rscript ${projectDir}/analysis/run.R config.auto.R
    echo "[RUN_ANALYSIS] done"
    """
}
