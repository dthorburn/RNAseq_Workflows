#!/bin/sh
#PBS -lwalltime=6:00:00
#PBS -lselect=1:ncpus=32:mem=32gb
#PBS -N X02_StringTie_Assembly
#PBS -j oe
#PBS -J 1-12

module load anaconda3/personal
source activate RNAseq

DATA_DIR=./02_Mapped_Reads
OUT_DIR=./03_Transcriptome

Annotations=./00_Reference_Index//Ref_Genome.gtf

Reads1=`sed ${PBS_ARRAY_INDEX}"q;d" /rds/general/user/dthorbur/home/tmstorage/ephemeral/sperm_agam_rna/00_R1_List.txt`
Newname=`echo $Reads1 | sed -e "s/_1.fastq//"`
Bam_File=`echo $Reads1 | sed -e "s/_1.fastq/.bam/"`

stringtie ${DATA_DIR}/${Bam_File} -p 30 -G ${Annotations} -o ${OUT_DIR}/${Newname}.gtf -l ${Newname} &&  echo "Finished assembling ${Newname} gtf: `date`"
