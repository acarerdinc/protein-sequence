from readdata import read_sequences, extract_kmers
import pyclust as pyclust
from Bio.SubsMat import MatrixInfo
import Bio.pairwise2 as pairwise2
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt
import numpy as np
import collections
import operator

# read protein sequences in a dict
sequences = read_sequences()

# extract kmers
sampling_size = 500
np.random.seed(42)
all_kmers = extract_kmers(sequences)
sub_kmers = np.random.choice(all_kmers, sampling_size)
print len(sub_kmers)

# compute the alignment score matrix
pp = 20
eps = np.finfo(np.float32).eps
kmer_alignments = np.zeros((len(sub_kmers), len(sub_kmers)))
for i in range(0, len(sub_kmers)):
    if pp != 0 and (i + 1) % (np.int(len(sub_kmers) / (100 / pp))) == 0:
        print('%{} completed.'.format(pp * (i + 1) /
                                      (np.int(len(sub_kmers) / (100 / pp)))))
    for j in range(0, len(sub_kmers)):
        alns = pairwise2.align.localds(sub_kmers[i], sub_kmers[j], MatrixInfo.blosum62,
                                       -50, -50)
        if len(alns) == 0:
            sim_score = 0
        else:
            sim_score = alns[0][2]
        kmer_alignments[i][j] = sim_score

# convert alignment scores to dissimilarities
# print kmer_alignments
kmer_diss = MinMaxScaler().fit_transform(kmer_alignments)
# print('after minmax norm')
# print kmer_diss
kmer_diss[kmer_diss == 0] += eps
kmer_diss = np.float32(kmer_diss)
# print('after eps')
# print kmer_diss
kmer_diss = -np.log(kmer_diss)
# print('after negative log')
# print kmer_diss

# bisecting k-means
bkm = pyclust.BisectKMeans(n_clusters=7, n_trials=20)
bkm.fit_predict(kmer_diss)
print('Bisecting Clustering Completed. ')

# clustering results
print collections.Counter(bkm.labels_)
bkm.tree_.show(line_type='ascii')
sorted_sse = sorted(bkm.sse_arr_.items(), key=operator.itemgetter(1))
print(sorted_sse)
# print(bkm.sse_arr_)
# print(bkm.sse_arr_.values())
sorted_sse_labels = [str(sorted_sse[i][0]) for i in range(0, len(sorted_sse))]
sorted_sse_y = [sorted_sse[i][1] for i in range(0, len(sorted_sse))]
plt.plot(sorted_sse_y)
# plt.xticks(sorted_sse_y, sorted_sse_labels, rotation='horizontal')
plt.show()
