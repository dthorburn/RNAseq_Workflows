# RNAseq Workflows:
# 1. STAR 2-pass
## Overview
This pipeline was developed to map RNAseq data using the [STAR](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/) (Spliced Transcripts Alignment to a Reference) 2-pass methodology. This pipeline is paired with variant calling using GATK in the repository [dthorburn/Genomic_Read_Processing](https://github.com/dthorburn/Genomic_Read_Processing). The pipeline was developed using [Nextfow](https://www.nextflow.io/) version 22.04.4 (version on Imperial HPC; date: 19/07/2022). 

**EDIT:** The Imperial HPC long nodes are no longer available. STAR is now more efficient to use as an array job. Please use the `starAlign.sh` script instead. The 15 samples listed in `00_File_List.txt` are a mix of paired and single end sequencing. You will need to change the values to define chromosomes if the number of samples change. 

### Usage

Update the paths in `00_File_List.txt`, update the working directory path in `starAlign.sh`, alter the `PBS_ARRAY_INDEX` statements to handle the chromosomes and samples correctly, and submit using `qsub starAlign.sh`

# 2. New Tuxedo Protocol (HISAT, StringTie, and Ballgown)
## Overview
This pipeline was developed to map RNAseq data using the "[new Tuxedo protocol](https://www.nature.com/articles/nprot.2016.095#Sec11)" (i.e., using [HISAT2](https://daehwankimlab.github.io/hisat2/), [StringTie](https://ccb.jhu.edu/software/stringtie/), and [Ballgown](https://www.bioconductor.org/packages/devel/bioc/vignettes/ballgown/inst/doc/ballgown.html)). This pipeline was developed using BASH, but on request I can update this to a nextflow workflow. 
