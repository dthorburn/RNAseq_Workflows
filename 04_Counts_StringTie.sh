#!/bin/sh
#PBS -lwalltime=6:00:00
#PBS -lselect=1:ncpus=32:mem=32gb
#PBS -N X04_StringTie_Counting
#PBS -j oe
#PBS -J 1-12

module load anaconda3/personal
source activate RNAseq



BAM_DIR=./02_Mapped_Reads
DATA_DIR=./03_Transcriptome
OUT_DIR=./04_Ballgown_Counts

Assembly=Assembly_ID
Annotations=./00_Reference_Index//Ref_Genome.gtf

Reads1=`sed ${PBS_ARRAY_INDEX}"q;d" /rds/general/user/dthorbur/home/tmstorage/ephemeral/sperm_agam_rna/00_R1_List.txt`
Newname=`echo $Reads1 | sed -e "s/_1.fastq//"`
Bam_File=`echo $Reads1 | sed -e "s/_1.fastq/.bam/"`

mkdir ${OUT_DIR}/${Newname}
stringtie -e -B -p 30 -G ${DATA_DIR}/${Assembly}.gtf -o ${OUT_DIR}/${Newname}/${Newname}_counts.gtf ${BAM_DIR}/${Bam_File} && echo "Finished counting ${Bam_File}: `date`"