# UGML
Analysis of NGS data from target gene exons capture sequencing    

UGML (Unveiling the Gene Mutational Landscape) is a Python script designed for Linux-based operating systems. The third-party programs required by UGML include fastp v0.23.2, flash v1.2.11, bowtie2, samtools, and Freebayes. UGML takes as input the randomly fragmented short-read sequencing results (in FASTQ format) from target gene exon capture regions and the target gene reference sequence (in FASTA format). Ultimately, it outputs the SNPs (.vcf) and amino acid substitution landscape (.txt) present in the target gene exons.
