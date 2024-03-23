import os
import sys

Codon_table =   {'TTT': 'F', 'CTT': 'L', 'ATT': 'I', 'GTT': 'V', 
                 'TTC': 'F', 'CTC': 'L', 'ATC': 'I', 'GTC': 'V', 
                 'TTA': 'L', 'CTA': 'L', 'ATA': 'I', 'GTA': 'V', 
                 'TTG': 'L', 'CTG': 'L', 'ATG': 'M', 'GTG': 'V', 
                 'TCT': 'S', 'CCT': 'P', 'ACT': 'T', 'GCT': 'A', 
                 'TCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A', 
                 'TCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A', 
                 'TCG': 'S', 'CCG': 'P', 'ACG': 'T', 'GCG': 'A', 
                 'TAT': 'Y', 'CAT': 'H', 'AAT': 'N', 'GAT': 'D', 
                 'TAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D', 
                 'TAA': '*', 'CAA': 'Q', 'AAA': 'K', 'GAA': 'E', 
                 'TAG': '*', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E', 
                 'TGT': 'C', 'CGT': 'R', 'AGT': 'S', 'GGT': 'G', 
                 'TGC': 'C', 'CGC': 'R', 'AGC': 'S', 'GGC': 'G', 
                 'TGA': '*', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G', 
                 'TGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'
                 }
class REF():
    def __init__(self,path):
        self.read = open(path,'r').read().split('\n')  # list.obj
        self.seq = (self.read)[1]        #string.obj
    def index(self):
        base_cnt = 0
        coding_base_cnt = 0
        base_index=[]
        AA_num=0
        exon_num=0
        switch = True
        for i in range(len(self.seq)):
            char = self.seq[i]
            base_cnt+=1    #the number of the base counting
            if char.islower():
                switch = True
                base_index.append([char,base_cnt,'Non-coding',0,0,'-'])
            elif char.isupper():
                if switch:
                    exon_num+=1
                switch = False
                coding_base_cnt+=1
                AA_num = (coding_base_cnt-1)//3+1
                codon_base_num = (coding_base_cnt-1)%3+1
                base_index.append([char,base_cnt,'coding',coding_base_cnt,AA_num,codon_base_num,f'ORF_{exon_num}'])
        return base_index
        
    def getcodon(self):
        string=''
        for i in range(len(self.seq)):
            char = self.seq[i]
            if char.isupper():
                string+=char
        codon_lst=[]
        cnt = 0
        i=0
        while i+3 < len(string):
            cnt+=1
            codon = string[i:i+3]
            AA = Codon_table[codon]
            codon_lst.append((AA,cnt,codon))
            i+=3
        return codon_lst   #[('M', 1, 'ATG'), ('S', 2, 'TCC'),..]
    
    def exon_index(self):
        bases_in_exon=[x for x in self.index() if 'coding' in x]
        codon_lst=self.getcodon()
        
        output=[]
        for element in bases_in_exon:
            dic={}
            AA_num=element[4]
            if (AA_num-1)>=len(codon_lst):
                break
            element[4]=codon_lst[AA_num-1]
            dic['REF_type']=element[0]
            dic['DNA_pos']=element[1]
            dic['CDS_pos']=element[3]
            dic['AA_info']=element[4]
            dic['CodonBaseNum']=element[5]
            dic['ORF']=element[6]
            output.append(dic)
        return output

class VCF():
    def __init__(self,path):
        self.init_index=[x for x in open(path,'r').read().split('\n') if (r'#' not in x and x != '')] #list.obj
    
    def preprocess(self):
        var_index=[]
        for i in self.init_index:
            dic={}
            info = i.split('\t')
            
            dic['REF_base'] = info[3]
            dic['DNA_pos_lst'] = [int(info[1])+x for x in range(len(info[3]))]
            dic['VAR_lst'] = info[4].split(',') 
            AO_buff = info[7].split(';')[5].split('=')[1]
            dic['DP']=int(info[7].split(';')[7].split('=')[1])
            dic['AO_lst']=[]
            if ',' in AO_buff:
                dic['AO_lst']=[int(x) for x in AO_buff.split(',')]
            else:
                dic['AO_lst'].append(int(AO_buff))
                    
            if len(dic['AO_lst'])!=0:
                if len(dic['REF_base'])==1:
                    dic['type']='SNP'
                    var_index.append(dic)
                elif len(dic['REF_base'])>1:
                    dic['type']='COMPLEX'
                    var_index.append(dic)
        return var_index
    
    def index(self):
        pre_idx=self.preprocess()
        output = []
        for x in pre_idx:
            for VAR,AO in zip(x['VAR_lst'],x['AO_lst']):
                dic={}
                dic['REF_type'] = x['REF_base']
                dic['DNA_pos_lst'] = x['DNA_pos_lst']
                dic['VAR_type'] = VAR
                dic['DP'] = x['DP']
                dic['AO'] = AO
                dic['type'] = x['type']
                output.append(dic)
        return output

def ExonVcfExtraction(Ref_path,Vcf_path):
    exon_idx=REF(Ref_path).exon_index()
    vcf_idx=VCF(Vcf_path).index()
    CDS_pos_range=[x['DNA_pos'] for x in exon_idx]
    output=[]
    for ele in vcf_idx:
        for pos in ele['DNA_pos_lst']:
            if pos in CDS_pos_range:
                ele['Base_info']=[x for x in exon_idx if pos ==x['DNA_pos']]
                output.append(ele)
                break
    return output       

def ORFSNPAnnotion(Ref_path,Vcf_path):
    exon_idx=REF(Ref_path).exon_index()
    vcf_idx=ExonVcfExtraction(Ref_path,Vcf_path)
    output=[]
    cnt=0
    for x in vcf_idx:
        cnt+=1
        if x['type'] == 'COMPLEX' and len(x['REF_type']) == len(x['VAR_type']):
            dic={}
            for REF_base,VAR_base,pos in zip(list(x['REF_type']),list(x['VAR_type']),x['DNA_pos_lst']):
                AA_buffer=[i for i in exon_idx if i['DNA_pos'] == pos]
                if len(AA_buffer) == 0:
                    continue
                CBN=AA_buffer[0]['CodonBaseNum']
                AA_info =AA_buffer[0]['AA_info']
                dic[pos]=(REF_base,VAR_base,AA_info,CBN)
            AA_lst=[]
            for info in dic.values():
                if info[2] not in AA_lst:
                    AA_lst.append(info[2])
            for i in AA_lst:
                VAR_buffer=list(i[2])  # ['A','T','G']
                for info in dic.values():
                    if i == info[2] and info[0]!=info[1]:
                        VAR_buffer[info[3]-1]=info[1]
                        
                VAR_codon = ''.join(VAR_buffer)
                VAR_AA = Codon_table[VAR_codon]
                if i[2] != VAR_codon:
                    if i[0] == VAR_AA:
                        output.append('\t'.join([i[0]+str(i[1])+VAR_AA,f"({str(i[2])}-->{VAR_codon})",str(x['AO'])+'/'+str(x['DP'])+f"({round(x['AO']/x['DP'],5)})",x['Base_info'][0]['ORF'],'Synonymous Mutations',f'Linkage_{cnt} Mutation']))        
                    elif i[0] != VAR_AA:
                        output.append('\t'.join([i[0]+str(i[1])+VAR_AA,f"({str(i[2])}-->{VAR_codon})",str(x['AO'])+'/'+str(x['DP'])+f"({round(x['AO']/x['DP'],5)})",x['Base_info'][0]['ORF'],'Non-synonymous Mutations',f'Linkage_{cnt} Mutation'])) 
                        
        elif x['type'] == 'SNP' and len(x['REF_type']) == len(x['VAR_type']):
            for REF_base,VAR_base,pos in zip(list(x['REF_type']),list(x['VAR_type']),x['DNA_pos_lst']):
                AA_buffer=[i for i in exon_idx if i['DNA_pos'] == pos]
                CBN=AA_buffer[0]['CodonBaseNum']
                
                REF_AA = AA_buffer[0]['AA_info'][0]
                AA_pos = AA_buffer[0]['AA_info'][1]
                REF_codon = AA_buffer[0]['AA_info'][2]
                
                VAR_codon_buffer = list(REF_codon)
                VAR_codon_buffer[CBN-1]=x['VAR_type']
                VAR_codon = ''.join(VAR_codon_buffer)
                VAR_AA =Codon_table[VAR_codon]
                if REF_AA == VAR_AA:
                    output.append('\t'.join([REF_AA+str(AA_pos)+VAR_AA,f"({REF_codon}-->{VAR_codon})",str(x['AO'])+'/'+str(x['DP'])+f"({round(x['AO']/x['DP'],5)})",x['Base_info'][0]['ORF'],'Synonymous Mutations',f'Linkage_{cnt}  Mutation']))
                elif REF_AA != VAR_AA:
                    output.append('\t'.join([REF_AA+str(AA_pos)+VAR_AA,f"({REF_codon}-->{VAR_codon})",str(x['AO'])+'/'+str(x['DP'])+f"({round(x['AO']/x['DP'],5)})",x['Base_info'][0]['ORF'],'Non-synonymous Mutations',f'Linkage_{cnt}  Mutation']))
    return output

def Make_table(REF_path,VCF_path,threashould=0):
    LST=[x.split('\t') for x in ORFSNPAnnotion(REF_path,VCF_path)]
    dic={}
    output=[]
    for i in range(len(LST)):
        if '\t'.join(LST[i][0:2]) not in dic.keys():
            dic['\t'.join(LST[i][0:2])] = [LST[i][2:5]]
        elif '\t'.join(LST[i][0:2]) in dic.keys():
            dic['\t'.join(LST[i][0:2])].append(LST[i][2:5])
    
    
    for x,y in dic.items():
        if len(y)>1:
            AO = 0
            for ele in y:
                AO+=int(ele[0].split('/')[0])
                DP = int(ele[0].split('/')[1].split('(')[0])
                ORF = ele[1]
                TYPE = ele[2]
            if AO/DP >= threashould:
                output.append('\t'.join([x,f"{str(AO)}/{str(DP)}({round(AO/DP,5)})",ORF,TYPE]))
        elif len(y) == 1:
            AO = 0
            for ele in y:
                AO+=int(ele[0].split('/')[0])
                DP = int(ele[0].split('/')[1].split('(')[0])
                ORF = ele[1]
                TYPE = ele[2]
            if AO/DP >= threashould:
                output.append('\t'.join([x,f"{str(AO)}/{str(DP)}({round(AO/DP,5)})",ORF,TYPE]))
    return output    

if __name__ == "__main__":

    MIN=0.01
    for i in range(len(sys.argv)):
        if '-R' in sys.argv[i]:
            ref=sys.argv[i+1]
        elif '-V' in sys.argv[i]:
            vcf=sys.argv[i+1]
        elif '-M' in sys.argv[i]:
            MIN=float(sys.argv[i+1])
        elif '-O' in sys.argv[i]:
            OUT=sys.argv[i+1]
    
    LST=Make_table(f"{ref}",
               f"{vcf}",threashould=MIN)
    Output=f'{OUT}'
    
    win=open(Output,'w',encoding='utf-8')
    print('\t'.join(['Mutations','Codon','Frequencies(%)','Location']),file=win)
    for x in LST:
        L=x.split('\t')
        L[2]=str(round(float(L[2].split('(')[1][:-1])*100,2))
        O='\t'.join(L[:5])
        print(O,file=win)
    win.close()