#!/usr/bin/env nextflow

/*
 * Pipeline developed for RNA alignment using STAR. 
 * Author: Miles Thorburn <d.thorburn@imperial.ac.uk>
 * Date last modified: 11/07/2022
 */

def helpMessage() {
  log.info """
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

"""
}

log.info """
==============================================================================================================================
                                        STAR Alignment Pipeline v1
==============================================================================================================================

Reference     : ${params.RefGen}
Annotations   : ${params.RefGTF}
Input         : ${params.InDir}
Aligned       : ${PWD}/02_Star_Pass2/

==============================================================================================================================
"""

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

if(!params.RefGen) {
  log.info"""
ERROR: No reference genome path provided! --RefGen /path/to/genome.fasta
==============================================================================================================================
  """
  helpMessage()
  exit 0
}

if(!params.InDir) {
  log.info"""
ERROR: No input directory provided! --InDir /path/to/reads/
==============================================================================================================================
  """
  helpMessage()
  exit 0
}

                                                            // =========================================================
                                                            // Setting the value channels (can be read unlimited times)
                                                            // =========================================================
ref_genome = file( params.RefGen, checkIfExists: true )
ref_annots = file( params.RefGTF, checkIfExists: true )
                                                            // =========================================================
                                                            // Step 1: Indexing Reference Genome
                                                            // =========================================================
// Setting up the chromosome channel
if( params.Chroms == "" ){
  // Defaulting to using all chromosomes
  chromosomes_ch = Channel
                      .from("AgamP4_2L", "AgamP4_2R", "AgamP4_3L", "AgamP4_3R", "AgamP4_X", "AgamP4_Y_unplaced", "AgamP4_UNKN")
  println "No chromosomes specified, using all major chromosomes: AgamP4_2L, AgamP4_2R, AgamP4_3L, AgamP4_3R, AgamP4_X, AgamP4_Y_unplaced, AgamP4_UNKN"
} else {
  // User option to choose which chromosome will be used
  // This worked with the following syntax nextflow run testing.nf --profile imperial --Chroms "AgamP4_3R,AgamP4_2L"
  chrs = params.Chroms.split(",")
  chromosomes_ch = Channel
                      .from(chrs)
  println "User defined chromosomes set: ${params.Chroms}"
}

if( params.Skip_STARPass1 == false ) {
  // Channels for reading in data
  if( params.Mode == "PE" ){
    Channel
      .fromFilePairs("${params.InDir}/*_{R1,R2,1,2}{.fastq.gz,.fq.gz,_001.fastq.gz,_001.fq.gz}")
      .ifEmpty { error "Cannot find any gzipped fastq files in ${params.InDir}" }
      .set { raw_fastqs }
  //    .map { ID, files -> tuple(ID, files[0], files[1])}
  } else if( params.Mode == "SE" ){
    Channel
      .fromPath("${params.InDir}/*{.fastq.gz,.fq.gz}")
      .ifEmpty { error "Cannot find any gzipped fastq files in ${params.InDir}" }
      .map { file -> tuple(file.baseName.replaceAll(".fastq.gz|.fq.gz|.fastq|.fq", ""), file) }
      .set { raw_fastqs }
  } else {
    log.info"""
ERROR: Incorrect mode input. Use either PE or SE.
==============================================================================================================================
    """
    helpMessage()
    exit 0
  }

  process STARPass1 {
    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries 3
    maxForks params.SP1_Forks

    tag { SampleID+"-"+chrom }

    publishDir(
      path: "${params.SP1Dir}",
      mode: 'copy',
      pattern: '*.tab'
    )

    executor = 'pbspro'
    clusterOptions = "-lselect=1:ncpus=${params.SP1_threads}:mem=${params.SP1_memory}gb -lwalltime=${params.SP1_walltime}:00:00"

    input:
    each chrom from chromosomes_ch
    set SampleID, file(reads) from raw_fastqs

    path ref_genome
    path ref_annots

    output:
    tuple SampleID, chrom, path("${SampleID}-${chrom}.SJ.out.tab"), path("reads*.fq.gz") into pass1_ch

    beforeScript 'module load star/2.7.1a; module load samtools/1.3.1'
    
    script:
    """
    mkdir -p STARref_${chrom}
    samtools faidx ${ref_genome}
    samtools faidx ${ref_genome} ${chrom} > ${chrom}.fasta
    samtools faidx ${chrom}.fasta && echo "Completed creating ${chrom} reference"

    STAR --runThreadN ${params.SP1_threads} \\
         --runMode genomeGenerate \\
         --genomeSAindexNbases ${params.SP1_IndexNbases} \\
         --genomeDir STARref_${chrom} \\
         --genomeFastaFiles ${chrom}.fasta ${params.SP1_ref_args} && echo "Completed generating reference index"

    STAR --runThreadN ${params.SP1_threads} \\
       --genomeDir STARref_${chrom} \\
       --readFilesCommand zcat \\
       --sjdbGTFfile ${ref_annots} \\
       --sjdbOverhang ${params.SP1_Overhang} \\
       --outFileNamePrefix ${SampleID}-${chrom}. \\
       --readFilesIn ${reads}  ${params.SP1_aln_args} && echo "Completed alingment pass 1"

    ## Cant save input symlinked files as output. So I'm copying them into the working directory to put into an output stream. This is ugly!
    COUNTER=0
    for file in ./*.gz
    do
      COUNTER=\$((COUNTER + 1))
      cp \${file} reads\${COUNTER}.fq.gz && echo "Copied reads\${COUNTER}.fq.gz"
    done
  """
  }
}
if( params.Skip_STARPass2 == false ) {
  if( params.Skip_STARPass1 ) {
    Channel
      .fromPath("${params.SP1Dir}/*.SJ.out.tab")
      .ifEmpty { error "Cannot find any index files in ${params.SP1Dir}" }
      .map { file -> tuple(file.baseName.replaceAll("-.*", ""), file.baseName.replaceAll("^.*-", "").replaceAll(".SJ.*", ""), file)}
      .set { temp_pass1_ch }

    if( params.Mode == "PE" ){
      Channel
        .fromFilePairs("${params.InDir}/*_{R1,R2,1,2}{.fastq.gz,.fq.gz,_001.fastq.gz,_001.fq.gz}")
        .ifEmpty { error "Cannot find any fastq files in ${params.InDir}" }
        .set { fastqs }
    //    .map { ID, files -> tuple(ID, files[0], files[1])}

    } else if( params.Mode == "SE" ){
      Channel
        .fromPath("${params.InDir}/*{.fastq.gz,.fq.gz,.fastq,.fq}")
        .ifEmpty { error "Cannot find any fastq files in ${params.InDir}" }
        .map { file -> tuple(file.baseName.replaceAll(".fastq.gz|.fq.gz|.fastq|.fq", ""), file) }
        .set { fastqs }
    } else {

    log.info"""
  ERROR: Incorrect mode input. Use either PE or SE.
  ==============================================================================================================================
      """
      helpMessage()
      exit 0
    }

    temp_pass1_ch
      .combine( fastqs, by: 0 )
      .map { sampleID, chrom, tab, reads -> tuple(sampleID, chrom, reads, tab) }
      .set { pass1_ch }
  }

  process STARPass2 {
    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries 3
    maxForks params.SP2_Forks

    tag { SampleID+"-"+chrom }

    publishDir(
      path: "${params.SP2Dir}",
      mode: 'copy',
    )

    executor = 'pbspro'
    clusterOptions = "-lselect=1:ncpus=${params.SP2_threads}:mem=${params.SP2_memory}gb -lwalltime=${params.SP2_walltime}:00:00"

    input:
    set SampleID, chrom, path(SJ_Tab), file(reads) from pass1_ch
    path ref_genome
    path ref_annots

    output:
    path("*.bam")
    path("*.bai")

    beforeScript 'module load star/2.7.1a; module load samtools/1.3.1'
    
    script:
    """
    mkdir -p STARref_${chrom}
    samtools faidx ${ref_genome}
    samtools faidx ${ref_genome} ${chrom} > ${chrom}.fasta
    samtools faidx ${chrom}.fasta && echo "Completed creating ${chrom} reference"

    STAR --runThreadN ${params.SP2_threads} \\
       --runMode genomeGenerate \\
       --genomeDir STARref_${chrom} \\
       --genomeFastaFiles ${chrom}.fasta \\
       --genomeSAindexNbases 11 \\
       --sjdbFileChrStartEnd ${SJ_Tab} ${params.SP2_ref_args} && echo "Completed generating reference index"

    STAR --runThreadN ${params.SP2_threads} \\
       --genomeDir STARref_${chrom} \\
       --readFilesCommand zcat \\
       --outFileNamePrefix ${SampleID}.${chrom}.2pass. \\
       --sjdbGTFfile ${ref_annots} \\
       --readFilesIn ${reads} \\
       --outSAMtype BAM SortedByCoordinate ${params.SP2_aln_args} && echo "Completed alingment pass 2"

    samtools index *.bam
    """
  }
}