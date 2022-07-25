#!/bin/sh
#PBS -lwalltime=72:00:00
#PBS -lselect=1:ncpus=64:mem=124gb
#PBS -N StarAlign_Array
#PBS -j oe
#PBS -J 1-75

module load star/2.7.1a
module load samtools/1.2

##############################################################
##                                                          ##
##                        Paramaters                        ##
##                                                          ##
##############################################################
PBS_ARRAY_INDEX=51 ##3,5,6,21,33,36,51
NSlots=62
Project_Dir=/rds/general/user/dthorbur/home/tmstorage/ephemeral/starAlign_PBSArray
ref_genome=/rds/general/project/tmstorage/live/agamp4_ref_genome/VectorBase-54_AgambiaePEST_Genome.fasta
ref_annots=/rds/general/project/tmstorage/live/agamp4_ref_genome/VectorBase-54_AgambiaePEST.gtf

##############################################################
##                                                          ##
##            Chromosome and Sample Distribution            ##
##                                                          ##
##############################################################

if [ ${PBS_ARRAY_INDEX} -le 15 ]
then
  Chrom="AgamP4_2R"
  Sample_INDEX=${PBS_ARRAY_INDEX}
elif [ ${PBS_ARRAY_INDEX} -le 30 ]
then
  Chrom="AgamP4_2L"
  let Sample_INDEX=${PBS_ARRAY_INDEX}-15
elif [ ${PBS_ARRAY_INDEX} -le 45 ]
then
  Chrom="AgamP4_3R"
  let Sample_INDEX=${PBS_ARRAY_INDEX}-30
elif [ ${PBS_ARRAY_INDEX} -le 60 ]
then
  Chrom="AgamP4_3L"
  let Sample_INDEX=${PBS_ARRAY_INDEX}-45
elif [ ${PBS_ARRAY_INDEX} -le 75 ]
then
  Chrom="AgamP4_X"
  let Sample_INDEX=${PBS_ARRAY_INDEX}-60
elif [ ${PBS_ARRAY_INDEX} -ge 76 ]
then
  echo "Too many samples"
  exit 0
fi

## Just for the test. 
## Chrom="AgamP4_Mt"

Forward_File=`sed ${Sample_INDEX}"q;d" ${Project_Dir}/00_File_List.txt`
SampleID=`basename ${Forward_File} | sed -e "s/_R1_001//" | sed -e "s/_1//" | sed -e "s/.fastq.gz//"`

## Handling the different file types. 
if echo $SampleID | grep -q "SRR7"
then
  Reads=${Forward_File}
else
  Reverse_File=`echo ${Forward_File} | sed -e "s/_R1_/_R2_/" | sed -e "s/_1/_2/"`
  Reads="${Forward_File} ${Reverse_File}"
fi

##############################################################
##                                                          ##
##                         Script                           ##
##                                                          ##
##############################################################

## Clunky as hell, but it works in an array like this. 
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "Starting analysis of ${SampleID} on ${Chrom}; `date`"
echo "Input Reads: ${Reads}"
echo "Output Folder: ${Project_Dir}/02_Output_Bams"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

mkdir -p ${Project_Dir}/01_Working_Dir/STARref_${Chrom}_${SampleID}_Pass1
cd ${Project_Dir}/01_Working_Dir/STARref_${Chrom}_${SampleID}_Pass1
samtools faidx ${ref_genome} ${Chrom} > ${Project_Dir}/01_Working_Dir/STARref_${Chrom}_${SampleID}_Pass1/${Chrom}.fasta
samtools faidx ${Project_Dir}/01_Working_Dir/STARref_${Chrom}_${SampleID}_Pass1/${Chrom}.fasta && echo "Completed creating ${Chrom} reference"

STAR --runThreadN ${NSlots} \
 --runMode genomeGenerate \
 --genomeSAindexNbases 11 \
 --genomeDir ${Project_Dir}/01_Working_Dir/STARref_${Chrom}_${SampleID}_Pass1 \
 --genomeFastaFiles ${Chrom}.fasta && echo "Completed generating reference index"

STAR --runThreadN ${NSlots} \
 --genomeDir ${Project_Dir}/01_Working_Dir/STARref_${Chrom}_${SampleID}_Pass1 \
 --readFilesCommand zcat \
 --sjdbGTFfile ${ref_annots} \
 --sjdbOverhang 149 \
 --outFileNamePrefix ${SampleID}_${Chrom}_1pass. \
 --readFilesIn ${Reads}  && echo "Completed alignment pass 1"

mkdir -p ${Project_Dir}/01_Working_Dir/STARref_${Chrom}_${SampleID}_Pass2
cd ${Project_Dir}/01_Working_Dir/STARref_${Chrom}_${SampleID}_Pass2
samtools faidx ${ref_genome} ${Chrom} > ${Project_Dir}/01_Working_Dir/STARref_${Chrom}_${SampleID}_Pass2/${Chrom}.fasta
samtools faidx ${Project_Dir}/01_Working_Dir/STARref_${Chrom}_${SampleID}_Pass2/${Chrom}.fasta 

STAR --runThreadN ${NSlots} \
 --runMode genomeGenerate \
 --genomeDir ${Project_Dir}/01_Working_Dir/STARref_${Chrom}_${SampleID}_Pass2 \
 --genomeFastaFiles ${Chrom}.fasta \
 --genomeSAindexNbases 11 \
 --sjdbFileChrStartEnd ${Project_Dir}/01_Working_Dir/STARref_${Chrom}_${SampleID}_Pass1/${SampleID}_${Chrom}_1pass.SJ.out.tab && echo "Completed generating reference index"

STAR --runThreadN ${NSlots} \
  --genomeDir ${Project_Dir}/01_Working_Dir/STARref_${Chrom}_${SampleID}_Pass2 \
  --readFilesCommand zcat \
  --outFileNamePrefix ${SampleID}_${Chrom}_2pass. \
  --sjdbGTFfile ${ref_annots} \
  --readFilesIn ${Reads} \
  --outSAMtype BAM SortedByCoordinate && echo "Completed alignment pass 2"

samtools index *.bam
rsync *.ba[mi] ${Project_Dir}/02_Output_Bams && echo "Finished tranferring files"

if [ -f  "${Project_Dir}/02_Output_Bams/${SampleID}_${Chrom}_2pass.Aligned.sortedByCoord.out.bam" ]
then
  rm -r ${Project_Dir}/01_Working_Dir/STARref_${Chrom}_${SampleID}_Pass1
  rm -r ${Project_Dir}/01_Working_Dir/STARref_${Chrom}_${SampleID}_Pass2
fi
