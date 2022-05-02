# SNPdragon
A high quality SNP calling pipeline for small to large microbial genomics projects.

## Installation

### Singularity (preferred)
`bash build.sh v1.0`[^docker]

### Conda
`conda env create --file docker/env_snpdragon_general.yml`  
`conda env create --file docker/env_snpdragon_freebayes.yml`  


## Usage
```
nextflow run main.nf \
--ref path/to/reference.fasta \
--out path/to/outdir \
--reads "path/to/reads/*_R{1,2}.f*.gz" \
-params-file params.yaml

Commands (required):
  --ref 	Reference file in fasta format  
  --out 	Absolute path to output directory  
  --reads 	Path to paired fastq files, must be enclosed in double quotes, 
		denote paired file pattern in curly braces (e.g. \*_R{1,2}.f\*.gz)
  -params-file  Additional parameters for variant filtering and output (params.yaml) 

params.yaml arguments:
  remove_clusters	Use sliding window approach to remove SNP clusters (default: false)  
  remove_cliffs	  	Remove SNPs occuring in read depth cliffs (default: true)  
  exclude_ref	  	Exclude reference from final matrix and alignment files (default: true)  
  contig_coverage	Minimum % breadth of coverage to include sample in output (default: 70)  
  mask_file	  	Optional .gff file to mask regions in reference (default: null)  
```

*For very large projects (1000s of isolates) where storage space for intermediate files may become an issue, main_v-cleanup.nf removes intermediate files[^intermediate] when they are no longer required while maintaining appropriate Nextflow re-entrancy behaviour. Requires absolute file paths.*

## Output
- 1_map/
  - \*.bam	: bwa-mem 
- 2_cov/
  - \*.pileup	: samtools depth  
- 3_call/
  - \*.bcf	: raw freebayes  
- 4_vcf/
  - \*.vcf	: filtered freebayes  
- 5_pkl/
  - \*.pkl	: pickled python Sample object  
- 6_msa/
  - core_snp.fasta	: core snp alignment  
  - core_stats.tsv	: per sample core genome size and % coverage over reference  
  - full_aln.fasta	: pseudo-genome alignment, low coverage positions denoted with 'N'  
  - full_snp.fasta	: whole-genome snp alignment, low coverage positions denoted with 'N'  
  - snp_dist.csv	: pairwise core snp distances  
  - snp_matrix.csv	: pairwise snp position x sample matrix  

[^docker]:  Docker for Mac (and possibly other platforms) limits resources to 2G by default.
  Increase this to at least 8G to build docker image.
[^intermediate]: Intermediate files in 1_map, 2_cov and 3_call are removed.
