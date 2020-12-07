#!/usr/bin/env python3

import os, gzip, itertools
from collections import Counter
# this is code which will parse FASTA files
# define what a header looks like in FASTA format
def isheader(line):
    return line[0] == '>'

def aspairs(f):
    seq_id = ''
    sequence = ''
    for header,group in itertools.groupby(f, isheader):
        if header:
            line = next(group)
            seq_id = line[1:].split()[0]
        else:
            sequence = ''.join(line.strip() for line in group)
            yield seq_id, sequence

def returnSum(myDict):

    sum = 0
    for i in myDict:
        sum = sum + myDict[i]

    return sum
url1="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2/cds/Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.ASM694v2.cds.all.fa.gz"
url2="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/mycobacterium_tuberculosis_h37rv/cds/Mycobacterium_tuberculosis_h37rv.ASM19595v2.cds.all.fa.gz"
file1="Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.ASM694v2.cds.all.fa.gz"
file2="Mycobacterium_tuberculosis_h37rv.ASM19595v2.cds.all.fa.gz"


if not os.path.exists(file1):
    os.system("curl -O %s"%(url1))

if not os.path.exists(file2):
    os.system("curl -O %s"%(url2))

n_gene1 = 0 #number of gene in file 1
n_gene2 = 0 #number of gene in file 2
L_gene1 = 0 #length of gene in file 1
L_gene2 = 0 #length of gene in file 2
L_gene2 = 0 #length of gene in file 2
G_number1 = 0
C_number1 = 0
G_number2 = 0
C_number2 = 0

codon1 = {}
# creating a dictionary for codon.
for first in {'A','T','C','G'}:
    for second in {'A','T','C','G'}:
        for third in {'A','T','C','G'}:
            new_codon = first+second+third
            new_dict = {new_codon:0}
            codon1.update(new_dict)

codon2 = {}
# creating a dictionary for codon.
for first in {'A','T','C','G'}:
    for second in {'A','T','C','G'}:
        for third in {'A','T','C','G'}:
         new_codon = first+second+third
         new__dict = {new_codon:0}
         codon2.update(new_dict)

codonF = {}
# creating a dictionary for codon.
for first in {'A','T','C','G'}:
    for second in {'A','T','C','G'}:
        for third in {'A','T','C','G'}:
            new_codon = first+second+third
            new_dict = {new_codon:0}
            codonF.update(new_dict)

with gzip.open(file1,"rt") as fh:
    seqs1 = aspairs(fh)
    for seq1 in seqs1:
        seqname1  = seq1[0]
        seqstring1= seq1[1]
        n_gene1 = n_gene1 + 1
        L_gene1 = L_gene1 + len(seqstring1)
        base = Counter(seqstring1)
		  L_gene1 = L_gene1 + len(seqstring1)
        base = Counter(seqstring1)
        G_species1 = G_number1 + base {'G'}]
        C_species1 = C_number1 + base {'C'}
        for i in range(0,len(seqstring1)-1,3):
            if seqstring1[i:i+3] in codon1:
                codon1[seqstring1[i:i+3]] += 1
            else:
                codon1[seqstring1[i:i+3]] = 1

with gzip.open(file2,"rt") as fj:
    seqs2 = aspairs(fj)
    for seq2 in seqs2:
        seqname2  = seq2[0]
        seqstring2= seq2[1]
        n_gene2 = n_gene2 + 1
        L_gene2 = L_gene2 + len(seqstring2)
        base = Counter(seqstring2)
        G_species2 = G_number2 + base['G']
        C_species2 = C_number2 + base['C']
        for i in range(0,len(seqstring2)-1,3):
            if seqstring2[i:i+3] in codon2:
                codon2[seqstring2[i:i+3]] += 1
            else:
                codon2[seqstring2[i:i+3]] = 1
       # print(seqname, " first 10 bases are ", seqstring[0:10])
for key in codonF:
    codonF[key] = str(round(codon1[key]/returnSum(codon1),4)) +  str(round(codon2[key]/returnSum(codon2),4))

print("The number of genes in species 1 is:",n_gene1)
print("The Length of genes in Species 1 is:",L_gene1)
print("The number of genes in species 2 is:",n_gene2)
print("The Length of genes in Species 2 is:",L_gene2)

print("The G+C percentage of species 1 is:",(G_species1 + C_species1)/L_gene1)
print("The G+C percentage of species 2 is:",(G_species2 + C_species2)/L_gene2)
print("The G+C percentage of total dataset is:",(G_species2 + C_species2 + G_species1 + C_species1 )/(L_gene2 + L_gene1))
print("The Codon numbers in file 1 are ",returnSum(codon1))
print("The Codon numbers in file 2 are ",returnSum(codon2))
print("The Codon frquency is", str(codonF))


