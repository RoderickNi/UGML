# UGML
Analysis of NGS data from target gene exons capture sequencing    

UGML (Unveiling the Gene Mutational Landscape) is a Python script designed for Linux-based operating systems. The third-party programs required by UGML include fastp v0.23.2, flash v1.2.11, bowtie2, samtools, and Freebayes. UGML takes as input the randomly fragmented short-read sequencing results (in FASTQ format) from target gene exon capture regions and the target gene reference sequence (in FASTA format). Ultimately, it outputs the SNPs (.vcf) and amino acid substitution landscape (.txt) present in the target gene exons.

## Dependencies
fastp v0.23.2 (https://github.com/OpenGene/fastp)    
FLASH v1.2.11 (https://github.com/dstreett/FLASH2)    
bowtie2 v2.3.5.1 (https://github.com/BenLangmead/bowtie2)    
samtools v1.3.1 (https://github.com/samtools/samtools)    
freebayes v1.3.7 (https://github.com/freebayes/freebayes)  

- Program execution flow   
![image execution flow](https://github.com/RoderickNi/UGML/blob/main/UGML_Program_execution_flow.png)

## Installation
- Conda environment    
```
conda create -n ugml python=3.10  # python version 3.6+
conda acitvate ugml
conda install fastp flash samtools freebayes bowtie2
```
- Get UGML
```
git clone https://github.com/RoderickNi/UGML.git
```

## Usage
```
python UGML.py --fq1  {RawReads_1.fq path}
               --fq2  {RawReads_2.fq path}
               --ref  {Reference.fa path}
               --od   {Output Directory path}
               --RT   {Read Type:  150 (PE150) or 250 (PE250) (default: 150)}
               --minF {Thresholds for SNP detection: minimum frequency (default: 0.01)}
               --minN {Thresholds for SNP detection: minimum number (default: 10)}
               --CPU  {CPU number for calculation (default: 3)}
```

## Example
```
python UGML.py --fq1 Raw_1.fq --fq2 Raw_2.fq --ref reference.fa --od OutputDir --RT 150 --minF 0.01 --minN 10 --CPU 3
```


