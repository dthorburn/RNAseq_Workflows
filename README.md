# RNAseq Workflows:
# 1. STAR 2-pass
## Overview
This pipeline was developed to map RNAseq data using the [STAR](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/) (Spliced Transcripts Alignment to a Reference) 2-pass methodology. This pipeline is paired with variant calling using GATK in the repository [dthorburn/Genomic_Read_Processing](https://github.com/dthorburn/Genomic_Read_Processing). The pipeline was developed using [Nextfow](https://www.nextflow.io/) version 22.04.4 (version on Imperial HPC; date: 19/07/2022). 

**EDIT:** The Imperial HPC long nodes are no longer available. STAR is now more efficient to use as an array job. Please use the `starAlign.sh` script instead. The 15 samples listed in `00_File_List.txt` are a mix of paired and single end sequencing. You will need to change the values to define chromosomes if the number of samples change. 

### Usage

Update the paths in `00_File_List.txt`, update the working directory path in `starAlign.sh`, alter the `PBS_ARRAY_INDEX` statements to handle the chromosomes and samples correctly, and submit using `qsub starAlign.sh`

# 2. New Tuxedo Protocol
## Overview
This pipeline was developed to map RNAseq data using the "[new Tuxedo protocol](https://www.nature.com/articles/nprot.2016.095#Sec11)" (i.e., using [HISAT2](https://daehwankimlab.github.io/hisat2/), [StringTie](https://ccb.jhu.edu/software/stringtie/), and [Ballgown](https://www.bioconductor.org/packages/devel/bioc/vignettes/ballgown/inst/doc/ballgown.html)). This pipeline was developed using BASH since for a high number of small jobs, a BASH array is more efficient at processing the samples than a nextflow script without the use of the deprecated Imperial HPC long nodes. On request I can update this to a nextflow workflow.

### Usage

First step is to create the `RNAseq` conda environment: `conda env create -f /path/to/RNAseq.yml`

Update the paths in `00_R1_List.txt` and update the `Reference`, `Annotations`, and `Assembly` variables to reflect the reference assembly in use. Next, update the number of jobs in the array `#PBS -J 1-X` where `X` is the number of samples in use in script numbers `01`, `02`, and `04`. Then submit each of the jobs in order from `00` to `04` using `qsub`. The job paramaters should be sufficient for most use cases. 
