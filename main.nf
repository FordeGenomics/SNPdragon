#!/usr/bin/env nextflow

/*
nextflow run main.nf \
--ref reference.fasta \
--out outdir \
--reads "reads/*_R{1,2}.f*.gz" \
-params-file params.yaml -profile test \
-with-timeline timeline.html \
-with-trace \
-with-report report.tsv
*/
nextflow.enable.dsl=1

log.info """\
reference       : $params.ref
out             : $params.out
reads           : $params.reads
remove_clusters : $params.remove_clusters
remove_cliffs   : $params.remove_cliffs
exclude_ref     : $params.exclude_ref
contig_coverage : $params.contig_coverage
mask_file       : $params.mask_file
"""

/* file object for reference and sample files */
ref_file = file(params.ref)
// samples_file = file(params.samples)

/*
Build genome index required for mapping
*/

process buildBwaIndex {
  publishDir "${params.out}/1_map", mode:'copy', overwrite:true
  tag "$ref_file.baseName"

  input:
  file ref from ref_file

  output:
  file "${ref}*" into mapping_ref_ch, freebayes_ref_ch, coverage_ref_ch

  script:
  """
  set +eu
  source activate snpdragon_general
  set -eu

  bwa index ${ref}
  samtools faidx ${ref}
  """
}

/*
Read in sample filenames
*/
// Channel
//     .fromPath("${params.samples}")
//     .splitText()
//     .map{ file(it.trim()).simpleName }
//     .set{ samples_ch }

/* Read in fastq files */

Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairs_ch }

// process trimming {
//   tag "${sample}"
//   publishDir "${params.out}/reads"
//   echo=true
//
//   input:
//   val sample from samples_ch
//
//   output:
//   set sample, file("${sample}_trimmed_R1.fastq.gz"), file("${sample}_trimmed_R2.fastq.gz") into mapping_ch
//
//   script:
//   """
//   set +eu
//   source activate snpdragon_general
//   set -eu
//
//   trimmomatic PE -threads 1 ${params.in}/${sample}_R1.fastq.gz ${params.in}/${sample}_R2.fastq.gz \
//   ${sample}_trimmed_R1.fastq.gz ${sample}_unpaired_R1.fastq.gz \
//   ${sample}_trimmed_R2.fastq.gz ${sample}_unpaired_R2.fastq.gz \
//   ILLUMINACLIP:/home.roaming/s4097594/db/adapters/trimmomatic.fa:1:30:11 LEADING:10 TRAILING:10 MINLEN:50 HEADCROP:10
//   """
// }

process mapping {
    label "multi"
    tag "${sample}"
    publishDir "${params.out}/1_map", mode:'copy', overwrite:true

    input:
    file ref from ref_file
    file ref_indexes from mapping_ref_ch
    tuple val(sample), file(reads) from read_pairs_ch

    output:
    tuple sample, file("${sample}.bam*") into freebayes_ch
    tuple sample, file("${sample}.bam*") into coverage_ch
    file("${sample}.bam*")

    script:
    """
    set +eu
    source activate snpdragon_general
    set -eu

    bwa mem \
    -R '@RG\\tID:1\\tSM:${sample}' \
    -t ${task.cpus} ${ref} \
    ${reads[0]} ${reads[1]} \
    | samtools sort -@ ${task.cpus} -O BAM \
    | samtools rmdup - ${sample}.bam

    samtools index ${sample}.bam

    """
    // bwa mem \
    // -R '@RG\\tID:1\\tSM:${sample}' \
    // -t ${task.cpus} ${ref} \
    // ${reads[0]} ${reads[1]} \
    // | samclip --max 10 --ref ${ref} \
    // | samtools sort -@ ${task.cpus} -O BAM \
    // | samtools rmdup - ${sample}.bam
}

process coverage {
    label "multi"
    tag "${sample}"
    publishDir "${params.out}/2_cov", mode:'copy', overwrite:true

    input:
    file ref from ref_file
    file ref_indexes from coverage_ref_ch
    tuple val(sample), file(bam) from coverage_ch

    output:
    tuple sample, file("${sample}.pileup") into pileup_ch
    file("${sample}.pileup")

    script:
    """
    set +eu
    source activate snpdragon_general
    set -eu

    parallel -k -j ${task.cpus} samtools depth -aa -q 10 -Q 60 -r {1} ${sample}.bam :::: <(fasta_generate_samtools_regions.py ${ref}.fai 100000 ) > ${sample}.pileup

    """
}

process freebayes {
    tag "${sample}"
    publishDir "${params.out}/3_call", mode:'copy', overwrite:true

    input:
    file ref from ref_file
    file ref_indexes from freebayes_ref_ch
    tuple val(sample), file(bam) from freebayes_ch

    output:
    tuple sample, file("${sample}.bcf*") into bcf_ch
    file("${sample}.bcf*")

    script:
    """
    set +eu
    source activate snpdragon_freebayes
    set -eu

    freebayes \
    -f ${ref} \
    --ploidy 2 \
    -P 0 \
    -F 0.1 \
    -q 10 \
    -m 10 \
    --min-alternate-count 2 \
    --min-repeat-entropy 1 \
    --genotype-qualities \
    ${sample}.bam \
    | vcfallelicprimitives -kg | bcftools view -O b -o ${sample}.bcf


    bcftools index ${sample}.bcf
    """
}

process samples {
    tag "${sample}"
    publishDir "${params.out}/4_vcf", pattern:'*vcf', mode:'copy', overwrite:true
    publishDir "${params.out}/5_pkl", pattern:'*pkl', mode:'copy', overwrite:true

    input:
    tuple val(sample), file(bcf), file(pileup) from bcf_ch.join(pileup_ch)

    output:
    file("${sample}.pkl") into pkl_ch
    val(sample) into process_ch
    val(sample) into sample_ch
    file("${sample}.pkl")
    file("*vcf")

    script:
    """
    set +eu
    source activate snpdragon_general
    set -eu

    pickle_samples.py -s ${sample} -c ${sample}.pileup -b ${sample}.bcf --remove_clusters ${params.remove_clusters} --remove_cliffs ${params.remove_cliffs}
    """
}

process process_results {
    label "big_memory"
    tag "${sample}"
    publishDir "${params.out}/6_msa/${sample}", mode:'copy', overwrite:true
    memory { 3.GB * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 5

    maxForks 10

    input:
    file ref from ref_file
    val samples from process_ch.collect()
    file pickles from pkl_ch.collect()
    val sample from sample_ch

    output:
    val(sample) into msa_ch
    file("*fasta") optional true into fasta_ch
    file("*core_stats.tsv") optional true into stats_ch
    file("*_snpmatrix.csv") optional true into mtx_ch
    file("*excluded.tsv") optional true into exclude_ch
    file("*fasta") optional true
    file("*csv") optional true
    file("*tsv") optional true

    """
    set +eu
    source activate snpdragon_general
    set -eu

    process_pickles.py -f "${ref}" -s "${sample}" -l "${samples}" -e ${params.exclude_ref} -c ${params.contig_coverage} --remove_clusters ${params.remove_clusters} --remove_cliffs ${params.remove_cliffs} 
    """
}

process concat {
  tag "Concatenating"
  publishDir "${params.out}/6_msa", mode:'copy', overwrite:true

  input:
  file ref from ref_file
  val samples from msa_ch.collect()
  file exclude_files from exclude_ch.collect().ifEmpty("No excluded samples")
  file fasta_files from fasta_ch.collect().ifEmpty("No included samples")
  file stat_files from stats_ch.collect().ifEmpty("No core stats")
  file mtx_files from mtx_ch.collect().ifEmpty("No SNP matrix")



  output:
  file("core_snp.fasta") optional true
  file("full_snp.fasta") optional true
  file("full_aln.fasta") optional true
  file("core_stats.tsv") optional true
  file("snp_dist.csv") optional true
  file("snp_matrix.csv") optional true
  file("excluded.tsv") optional true

  """
  set +eu
  source activate snpdragon_general
  set -eu

  concat_msa.py -f "${ref}" -l "${samples}" -o ${params.out}
  """
}
