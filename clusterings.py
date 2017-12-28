import sklearn.cluster as cluster
from Bio import SeqIO, pairwise2
import re
from Bio.SubsMat import MatrixInfo
import Bio.pairwise2 as pairwise2
from Bio.pairwise2 import format_alignment
import numpy as np

fh = open("uniprot-part.fasta.txt", 'r')
line = fh.readline()
meta = ''
sequence = ''
metalist = []
seqlist = []
idlist = []
n_samples = 5

while line:
    line = line.rstrip('\n')
    # print(line)
    if '>' in line:
        if sequence.__len__() > 1:
            seqlist.append(sequence)
        sequence = ''
        meta = line
        metalist.append(meta)
        seqobj = re.search(r'>sp.(\w*)', meta)
        idlist.append(seqobj.group(1))
    else:
        sequence = sequence + line
    line = fh.readline()
seqlist.append(sequence)

# print(metalist)
# print(seqlist)
# print(idlist)
# print(seqlist)


BLOSUM62_MTRX = MatrixInfo.blosum62


abbrs = {'ALA': 'A',
         'ARG': 'R',
         'ASN': 'N',
         'ASP': 'D',
         'CYS': 'C',
         'GLU': 'E',
         'GLN': 'Q',
         'GLY': 'G',
         'HIS': 'H',
         'ILE': 'I',
         'LEU': 'L',
         'LYS': 'K',
         'MET': 'M',
         'PHE': 'F',
         'PRO': 'P',
         'SER': 'S',
         'THR': 'T',
         'TRP': 'W',
         'TYR': 'Y',
         'VAL': 'V'}


from skbio import Protein
km = []
for i in range(len(seqlist)):
    s = Protein(seqlist[i])
    for kmer in s.iter_kmers(5, overlap=True):
        km.append(str(kmer))

# Subsample data
import random
random.seed(42)
print(km.__len__())
km = random.sample(km, n_samples)

print(km.__len__())

#km = ['AMNQHRV', 'VNKLMAV', 'GPARQHH']
num_seqs = (len(km))
sc = []
import numpy
mat = numpy.zeros(shape=(km.__len__(),km.__len__()))
print(mat)
#scores = [[0 for i in range(num_seqs)] for j in range(1, num_seqs)]
for i in range(num_seqs):
    for j in range(i, num_seqs):
        alns = pairwise2.align.localds(str(km[i]), str(km[j]), BLOSUM62_MTRX, -5, -1)   # TODO: sometimes returns empty, what does it mean?
        print('km[i]={}, km[j]={}, alns={}'.format(km[i], km[j], alns))
        #print (alns)
        #print([x[2] for x in alns][0])
#        scores.append([x[2] for x in alns][0])
 #       sc.append([x[2] for x in alns][0])
        if len(alns) == 0:
            mat[i][j] = 0
        else:
            mat[i][j]= ([x[2] for x in alns][0])
        #print(scores)
print('matrix:')
print(mat)
#print (sc)

mat = 1.0 / mat
mat[mat==np.Inf] = 0

print(mat)

kmeans = cluster.KMeans(2)
results = kmeans.fit(mat)
print(results)

labels = results.labels_
cls = [[] for i in range(2)]

for i in range(0, num_seqs):
    cls[labels[i]].append(km[i])
print(cls)

'''for i in range(0, len(cls)):
    output_handle = open("c." + str(i) + ".fasta", "w")
    SeqIO.write(cls[i], output_handle, "fasta")
    output_handle.close()'''