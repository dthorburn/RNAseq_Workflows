#!/bin/sh
#PBS -lwalltime=6:00:00
#PBS -lselect=1:ncpus=32:mem=32gb
#PBS -N X01_HISAT2_Mapping
#PBS -j oe
#PBS -J 1-12

module load hisat/2.0.4
module load python/2.7.3
module load samtools/1.3.1

DATA_DIR=./01_Raw_Data
OUT_DIR=./02_Mapped_Reads

Reference=./00_Reference_Index/Ref_Genome.fasta
Annotations=./00_Reference_Index//Ref_Genome.gtf
Ref_Index=./00_Reference_Index/genome_tran

Reads1=`sed ${PBS_ARRAY_INDEX}"q;d" /rds/general/user/dthorbur/home/tmstorage/ephemeral/sperm_acol_rna/00_R1_List.txt`
Reads2=`echo $Reads1 | sed -e "s/_1.fastq/_2.fastq/"`
Newname=`echo $Reads1 | sed -e "s/_1.fastq//"`
echo $Reads1 $Reads2 $Newname

hisat2 -p 30 --dta -x ${Ref_Index} -1 ${DATA_DIR}/${Reads1} -2 ${DATA_DIR}/${Reads2} -S ${OUT_DIR}/${Newname}.sam && echo "Finished mapping ${Newname}: `date`"

samtools sort --threads 30 -o ${OUT_DIR}/${Newname}.bam ${OUT_DIR}/${Newname}.sam && rm ${OUT_DIR}/${Newname}.sam 
samtools index ${OUT_DIR}/${Newname}.bam && echo "Finished indexing ${Newname} bam: `date`"