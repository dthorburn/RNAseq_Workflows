#!/bin/sh
#PBS -lwalltime=6:00:00
#PBS -lselect=1:ncpus=8:mem=32gb
#PBS -N X03_StringTie_Merge
#PBS -j oe

module load anaconda3/personal
source activate RNAseq

DATA_DIR=./03_Transcriptome

Assembly=Assembly_ID
Annotations=./00_Reference_Index//Ref_Genome.gtf

cd $DATA_DIR

stringtie --merge -p 8 -G ${Annotations} -o ${DATA_DIR}/${Assembly}.gtf $DATA_DIR/mergelist.txt && echo "Finished merging: `date`"
gffcompare -r ${Annotations} -G -o ${Assembly}_merged ${DATA_DIR}/${Assembly}.gtf && echo "Finished comparing: `date`"