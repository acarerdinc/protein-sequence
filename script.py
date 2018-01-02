import pyclust as pyclust
import scipy
import sklearn.cluster as cluster
import re
from Bio.SubsMat import MatrixInfo
import Bio.pairwise2 as pairwise2
import numpy as np
from skbio import Protein
import matplotlib.pyplot as plt


def km_mat(km):
    sc = []
    num_seqs = len(km)
    BLOSUM62_MTRX = MatrixInfo.blosum62
    mat = np.zeros(shape=(km.__len__(), km.__len__()))

    for i in range(num_seqs):
        for j in range(i, num_seqs):
            alns = pairwise2.align.localds(str(km[i]), str(km[j]), BLOSUM62_MTRX, -5, -1)
            print('km[i]={}, km[j]={}, alns={}'.format(km[i], km[j], alns))
            # print("i:",i)
            if len(alns) == 0:
                mat[i][j] = 0
            else:
                mat[i][j] = ([x[2] for x in alns][0])

    print('matrix:')
    print(mat)

    mat = 1.0 / mat
    mat[mat == np.Inf] = 0
    print(mat)
    return mat


sse_arr = []
steps = []

fh = open("C:\Python27\seq.fasta", 'r')
line = fh.readline()
meta = ''
sequence = ''
metalist = []
seqlist = []
idlist = []


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
print(seqlist)
# print(idlist)


km = []
for i in range(len(seqlist)):
    s = Protein(seqlist[i])
    for kmer in s.iter_kmers(5, overlap=True):
        km.append(str(kmer))
# print(km)
# import sys

# Subsample data
km = km[0:10]
num_seqs = len(km)

mat = km_mat(km)

for k in range(2, len(km)):
    kmeans = cluster.KMeans(k)
    results = kmeans.fit(mat)
    # print('result:\n', results)
    # print('SSE:\n', results.inertia_)
    # print(results.labels_)
    labels = results.labels_
    cls = [[] for i in range(k)]

    for i in range(0, num_seqs):
        cls[labels[i]].append(km[i])
    print(cls)
    sse_arr.append(results.inertia_)
    steps.append(k-1)
    print('k:')
    print(k)
	
	
	
print("---------------------------------------------------------------------------")


plt.plot(steps, sse_arr)
plt.xlabel("k", fontsize=14)
plt.ylabel("SSE", fontsize=14)
plt.show()


#Bisecting Kmeans

bkm = pyclust.BisectKMeans(n_clusters=3)
bkm.fit(mat)
labels = bkm.labels_
print(labels)

cls = [[] for i in range(k+10)]

for i in range(0, num_seqs):
    cls[labels[i]].append(km[i])
print(cls)


#SSE
# sum1 = 0
# sse1 = 0
# sum2 = 0
# sse2 = 0
#
# for i in range(len(cls[1])):
#     sum1 =+ cls[1][i]
# sum1 / int(len(cls[1]))
# for i in range(len(cls[1])):
#     sse1 =+ (cls[1][i]- sum)**2
#
# print(sse1)
#
#
# for i in range(len(cls[2])):
#     sum2 =+ cls[2][i]
# sum2 / int(len(cls[2]))
# for i in range(len(cls[2])):
#     sse2 =+ (cls[2][i]- sum)**2
# print(sse2)
#
#
# mat_bisect = numpy.zeros(shape=(cls[1].__len__(),cls3[1].__len__()))
# for i in range (len(cls[1])):
#     for j in range(i, len(cls[1])):
#         alns_bisect = pairwise2.align.localds(str(cls[1][i]), str(cls[1][j]), BLOSUM62_MTRX, -5, -1)
#         if len(alns_bisect) == 0:
#             mat[i][j] = 0
#         else:
#             mat_bisect[i][j] = ([x[2] for x in alns_bisect][0])
#             # print(scores)
# print('matrix:')
# print(mat_bisect)


# print("---------------------------------------------------------------------------")
