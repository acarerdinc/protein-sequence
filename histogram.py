import re
import numpy as np
import matplotlib.pyplot as plt
import pylab as plt

fh = open("uniprot-2006part.fasta.txt", 'r')
line = fh.readline()
meta = ''
sequence = ''
metalist = []
seqlist = []
idlist = []
n_samples = 10

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

#print(metalist)
#print(seqlist)
#print(idlist)
#print(seqlist)

frequency = []
for i in range (len(seqlist)):
    frequency.append(len(seqlist[i]))
#print(frequency)

id = []
for i in range(len(idlist)):
    id.append(i)
#print(id)
#print(len(id))

frequency.sort(reverse=True)
plt.bar(id, frequency, align='center')
plt.xticks([])
plt.xlabel('Protein ID', fontsize=16)
plt.ylabel('Length', fontsize=16)
#plt.xticks(id, idlist)
plt.show()
