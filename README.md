# RNAseq Read Mapping
## Overview
This pipeline was developed to map RNAseq data using the [STAR](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/) (Spliced Transcripts Alignment to a Reference) 2-pass methodology. This pipeline is paired with variant calling using GATK in the repository [dthorburn/Genomic_Read_Processing](https://github.com/dthorburn/Genomic_Read_Processing). The pipeline was developed using [Nextfow](https://www.nextflow.io/) version 22.04.4 (version on Imperial HPC; date: 19/07/2022). 

**EDIT:** The Imperial HPC long nodes are no longer available. STAR is now more efficient to use as an array job. Please use the `starAlign.sh` script instead. The 15 samples listed in `00_File_List.txt` are a mix of paired and single end sequencing. You will need to change the values to define chromosomes if the number of samples change. 

### Usage
Below is the help message from `STAR_Align.nf` including instructions on how to run the pipeline:
```
Usage:
  This pipelines was developed to map RNA reads to a reference genome using the STAR 2-pass methodology.

To use, there are 3 steps:
  1. Update project directory path in STAR_Align.sh
  2. Add required arguments listed below
  3. Submit pipeline coordinator using qsub STAR_Align.sh

If you require available HPC jobs for alternative scripts lower job concurrency options.

Required arguments:
  --RefGen                                    Path to reference fasta. Usage: '--RefGen /path/to/genome.fasta'
  --RefGTF                                      
  --InDir                                     Path to input gzipped fastq directory. Required even if skipping pass 1.
  --Mode                                      Paired-end or single-end input. Usage: '--Mode PE' or '--Mode SE'

Optional arguments:
  -w                                          Path to nextflow working directory. (Default: ./work)
  --help                                      Show this message
  --Chroms                                    User defined chromosome selection. (Default: all major LGs in AgamP4).
                                              Usage: '--Chroms "AgamP4_3R,AgamP4_X,AgamP4_Mt'. Selection much be comma
                                              delimited with no spaces and match the contig names in the fasta index.
  --SP1_ref_args                              Optional arguments for STAR genomeGenerate
  --SP1_aln_args                              Optional arguments for STAR mapping
  --SP2_ref_args                              Optional arguments for STAR genomeGenerate pass 2
  --SP2_aln_args                              Optional arguments for STAR mapping pass 2

Concurrency options:                          Imperial HPC permits a maximum of 50 jobs per user. 
  --SP1_Forks                                 Default: 15
  --SP2_Forks                                 Default: 15

Debugging arguments:
  --Skip_STARPass1                            Skip STAR genomeGenerate and mapping pass 1
  --Skip_STARPass2                            Skip STAR genomeGenerate and mapping pass 2
  --SP1_threads                               Number of threads for each subprocess - swap BP for process any acronym to
                                              alter other processes. (i.e., SP2_walltime = 24)
  --SP1_memory                                Number of Gb of memory for each subprocess
  --SP1_walltime                              Number of hours for each subprocess (72 is maximum)                            
```
