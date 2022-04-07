#!/usr/bin/env nextflow

/*
nextflow run main_v-cleanup.nf --ref /path/to/reference.fasta --out /path/to/outdir --reads "/path/to/reads/*_R{1,2}.f*.gz" -with-timeline timeline.html -with-trace -with-report report.tsv -profile slurm_cluster

*/

log.info """\
reference     : $params.ref
out           : $params.out
reads         : $params.reads
"""

/* file object for reference and sample files */
ref_file = file(params.ref)
// samples_file = file(params.samples)

/*
Build genome index required for mapping
*/

process buildBwaIndex {
  publishDir "${params.out}/1_map", overwrite:true
  tag "$ref_file.baseName"

  input:
  file ref from ref_file

  output:
  file "${ref}*" into mapping_ref_ch, freebayes_ref_ch, coverage_ref_ch

  script:
  """
  set +eu
  source activate snp_pipeline_general
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
//   source activate snp_pipeline_general
//   set -eu
//
//   mkdir -p ${params.out}/reads
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
    publishDir "${params.out}/1_map", pattern:'*bam*', overwrite:true

    input:
    file ref from ref_file
    file ref_indexes from mapping_ref_ch
    tuple val(sample), file(reads) from read_pairs_ch

    output:
    tuple sample, file("${sample}.mapping") into freebayes_ch
    tuple sample, file("${sample}.mapping") into coverage_ch
    file("${sample}.bam*")

    script:
    """
    set +eu
    source activate snp_pipeline_general
    set -eu

    mkdir -p ${params.out}/1_map

    echo Mapping ${sample} `date '+%d/%m/%Y_%H:%M:%S'`

    bwa mem \
    -R '@RG\\tID:1\\tSM:${sample}' \
    -t ${task.cpus} ${ref} \
    ${reads[0]} ${reads[1]} \
    | samtools sort -@ ${task.cpus} -O BAM \
    | samtools rmdup - ${sample}.bam

    samtools index ${sample}.bam

    touch ${sample}.mapping
    """
}

process coverage {
    label "multi"
    tag "${sample}"
    publishDir "${params.out}/2_cov", pattern:'*pileup', overwrite:true

    input:
    file ref from ref_file
    file ref_indexes from coverage_ref_ch
    tuple val(sample), file(mapping) from coverage_ch

    output:
    tuple sample, file("${sample}.coverage") into pileup_ch
    file("${sample}.pileup")

    script:
    """
    set +eu
    source activate snp_pipeline_general
    set -eu

    mkdir -p ${params.out}/2_cov

    parallel -k -j ${task.cpus} samtools depth -aa -q 10 -Q 60 -r {1} ${params.out}/1_map/${sample}.bam :::: <(/${params.script_dir}/fasta_generate_samtools_regions.py ${ref}.fai 100000 ) > ${sample}.pileup

    touch ${sample}.coverage

    """
    // samtools mpileup -f ${ref} -a -q 1 ${bam} -o ${sample}.pileup
}

process freebayes {
    tag "${sample}"
    publishDir "${params.out}/3_call", pattern:'*bcf*', mode:'copy', overwrite:true

    input:
    file ref from ref_file
    file ref_indexes from freebayes_ref_ch
    tuple val(sample), file(mapping) from freebayes_ch

    output:
    tuple sample, file("${sample}.freebayes") into bcf_ch
    file("${sample}.bcf*")

    script:
    """
    set +eu
    source activate snp_pipeline_freebayes
    set -eu

    mkdir -p ${params.out}/3_call

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
    ${params.out}/1_map/${sample}.bam \
    | vcfallelicprimitives -kg | bcftools view -O b -o ${sample}.bcf

    bcftools index ${sample}.bcf

    touch ${sample}.freebayes
    """
}

process samples {
    tag "${sample}"
    publishDir "${params.out}/4_vcf", pattern:'*vcf', mode:'copy', overwrite:true
    publishDir "${params.out}/5_pkl", pattern:'*pkl', mode:'copy', overwrite:true

    input:
    tuple val(sample), file(freebayes), file(coverage) from bcf_ch.join(pileup_ch)

    output:
    file("${sample}.pkl") into pkl_ch
    val(sample) into process_ch
    val(sample) into sample_ch
    val(sample) into cleanup_ch
    file("${sample}.pkl")
    file("*vcf")

    script:
    """
    set +eu
    source activate snp_pipeline_general
    set -eu

    mkdir -p ${params.out}/4_vcf
    mkdir -p ${params.out}/5_pkl

    pickle_samples.py -s ${sample} -c ${params.out}/2_cov/${sample}.pileup -b ${params.out}/3_call/${sample}.bcf -r

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
    file("*fasta") into fasta_ch
    file("*tsv") into stats_ch
    file("*_snpmatrix.csv") into mtx_ch
    val(sample) into msa_ch
    file("*fasta")
    file("*csv")

    """
    set +eu
    source activate snp_pipeline_general
    set -eu
    
    mkdir -p ${params.out}/6_msa/${sample}

    process_pickles.py -f "${ref}" -s "${sample}" -l "${samples}" -o ${params.out} -e -r
    """
}

process concat {
  tag "Concatenating"
  publishDir "${params.out}/6_msa", mode:'copy', overwrite:true

  input:
  file ref from ref_file
  file fasta_files from fasta_ch.collect()
  file stat_files from stats_ch.collect()
  file mtx_files from mtx_ch.collect()
  val samples from msa_ch.collect()

  output:
  file("core_snp.fasta")
  file("full_snp.fasta")
  file("full_aln.fasta")
  file("core_stats.tsv")
  file("snp_dist.csv")
  file("snp_matrix.csv")
  """
  set +eu
  source activate snp_pipeline_general
  set -eu

  concat_msa.py -f "${ref}" -l "${samples}" -o ${params.out}
  """
}

process cleanup {
  tag "${sample}"

  input:
  val sample from cleanup_ch

  """
  echo | tee `readlink ${params.out}/1_map/${sample}.bam*`
  echo | tee `readlink ${params.out}/2_cov/${sample}.pileup`
  echo | tee `readlink ${params.out}/3_call/${sample}.bcf*`
  """
}
