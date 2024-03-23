import os
import sys


def UGML(RawRead_1,RawRead_2,ref,OutDir,CPU_num,ReadType,MinFqcy,MinNumb):
    print()
    print("="*50)
    print("Step1: Quality Control")
    print("="*50)
    os.system(f"rm -rf {OutDir}")
    os.system(f"mkdir {OutDir}")
    C1=os.path.join(OutDir,'CleanRead1.fq')
    C2=os.path.join(OutDir,'CleanRead2.fq')
    os.system(f'''
        fastp --in1 {RawRead_1} --in2 {RawRead_2} --out1 {C1} --out2 {C2} --n_base_limit 0 --unqualified_percent_limit 0 --qualified_quality_phred 30 --length_required 100 --correction
        rm fastp*
    ''')
    
    print()
    print("="*50)
    print("Step2: Overlap Assembly")
    print("="*50)
    OVLP=os.path.join(OutDir,'overlap')
    os.system(f'''
        flash {C1} {C2} -d {OVLP} -t {CPU_num} -p 33 -m 10 -M {ReadType}
    ''')
    
    print()
    print("="*50)
    print("Step3: Map to Referenence")
    print("="*50)
    OVfasta=os.path.join(OVLP,'out.extendedFrags.fastq')
    os.system(f'''
        bowtie2-build {ref} Ref_index
        bowtie2 -p {CPU_num} -x Ref_index -U {OVfasta} | samtools sort -m 16G -o MapRst.bam
    ''')
    
    print()
    print("="*50)
    print("Step4: Extraction of mapped reads")
    print("="*50)
    os.system(f'''
        samtools view -bF 4 MapRst.bam > MappedReads.bam
        samtools faidx {ref}
        samtools index -b MappedReads.bam
    ''')
    print("Done!")
    
    print()
    print("="*50)
    print("Step5: SNP Calling")
    print("="*50)
    VAR=os.path.join(OutDir,'var.vcf')
    os.system(f'''
       freebayes -f {ref} -F {MinFqcy} -C {MinNumb} --pooled-continuous MappedReads.bam > {VAR}
    ''')
    print("Done!")

    print()
    print("="*50)
    print("Step6: Getting Amino acid substitution landscape")
    print("="*50)
    OUTAAvar=os.path.join(OutDir,'AA_var.txt')
    os.system(f'''
       python SubLand.py -R {ref} -V {VAR} -O {OUTAAvar} -M {MinFqcy}
    ''')
    print("Done!")
    
    os.system(f'''
        rm *.bt2
        rm *.bam
        rm *.bai
    ''')
    


if __name__ == "__main__":

    
    CPU_num=3
    ReadType=150
    MinFqcy=0.01
    MinNumb=10
    for i in range(len(sys.argv)):
        if '--fq1' == sys.argv[i]:       # RawReads_1
            RawRead_1=sys.argv[i+1]
        elif '--fq2' == sys.argv[i]:     # RawReads_2
            RawRead_2=sys.argv[i+1]
        elif '--ref' == sys.argv[i]:     # Reference
            ref = sys.argv[i+1]
        elif '--od' == sys.argv[i]:      # Output Dir Path
            OutDir = sys.argv[i+1]
        elif '--CPU' == sys.argv[i]:     # CPU number for calculation
            CPU_num=int(sys.argv[i+1])
        elif '--RT' == sys.argv[i]:      # Read Type :  150 or 250
            ReadType=int(sys.argv[i+1])
        elif '--minF' == sys.argv[i]:    # Thresholds for SNP detection: minimum frequency
            MinFqcy = sys.argv[i+1]
        elif '--minN' == sys.argv[i]:    # Thresholds for SNP detection: minimum number
            MinNumb = sys.argv[i+1]

    
    UGML(RawRead_1,RawRead_2,ref,OutDir,CPU_num,ReadType,MinFqcy,MinNumb)

