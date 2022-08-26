#!/bin/sh
#PBS -lwalltime=6:00:00
#PBS -lselect=1:ncpus=32:mem=256gb
#PBS -N X00_HISAT2_Index
#PBS -j oe

mkdir 00_Reference_Index
mkdir 01_Raw_Data
mkdir 02_Mapped_Reads
mkdir 03_Transcriptome
mkdir 04_Ballgown_Counts

module load hisat/2.0.4
module load python/2.7.3
cd 00_Reference_Index

Reference=Ref_Genome.fasta
Annotations=Ref_Genome.gtf

echo "Starting: `date`"

hisat2_extract_splice_sites.py ${Annotations} > genome.ss && echo "Finished extracting splice sites"
hisat2_extract_exons.py ${Annotations} > genome.exon && echo "Finished extracting exons"

hisat2-build -p 30 --exon genome.exon --ss genome.ss ${Reference} genome_tran && echo "All done: `date`"