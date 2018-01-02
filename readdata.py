""" Author: Acar Erdinc
"""

from os import listdir
from os.path import isfile, join
import fnmatch
import re


def read_sequences(filename="Spmap/uniprot-all.fasta"):
    """
    Reads protein sequences and ids.
    :return: A dictionary of protein sequences. Protein IDs are keys.
    """
    fh = open(filename, 'r')

    sequences = {}
    sequence = ''
    line = fh.readline()
    while line:
        line = line.rstrip('\n')
        if '>' in line:
            if len(sequence) > 1:
                sequences[seq_id] = sequence
            sequence = ''
            meta = line
            seq_obj = re.search(r'>sp.(\w*)', meta)
            seq_id = seq_obj.group(1)
        else:
            sequence = sequence + line
        line = fh.readline()
    if len(sequence) > 1:
        sequences[seq_id] = sequence

    return sequences


def read_data(level=0, length_limit=None):
    """
    Reads sequence and hierarchical enzyme clustering labels.
    :param level: Level of the clustering from general to specific
    :return: (x_train, y_train, x_test, y_test). Training and test data and targets
    """
    sequences = read_sequences()
    x_train = []
    y_train = []
    x_test = []
    y_test = []
    duplicates = []

    path_positive_train = 'Spmap/positiveTrain'
    path_positive_test = 'Spmap/positiveTest'

    if level == 0:
        regex = fnmatch.translate('[0-9].-*.ids')
    elif level == 1:
        regex = fnmatch.translate('[0-9].[0-9].-*.ids')
    re_obj = re.compile(regex)

    # Read train data
    file_list = [f_name for f_name in listdir(path_positive_train)
                 if isfile(join(path_positive_train, f_name))]
    for f_name in file_list:
        if re_obj.match(f_name):
            y = f_name[0:(level * 2 + 1)]
            with open(join(path_positive_train, f_name)) as f:
                content = f.readlines()  # Read file content
            content = [x.strip() for x in content]
            for p in content:
                if p in sequences:
                    if sequences[p] in x_train:
                        duplicates.append(p)
                    if length_limit is None or len(sequences[p]) <= length_limit:
                        x_train.append(sequences[p])
                        y_train.append(y)

    # Read test data
    file_list = [f_name for f_name in listdir(path_positive_test)
                 if isfile(join(path_positive_test, f_name))]
    for f_name in file_list:
        if re_obj.match(f_name):
            y = f_name[0:(level * 2 + 1)]
            with open(join(path_positive_test, f_name)) as f:
                content = f.readlines()  # Read file content
            content = [x.strip() for x in content]
            for p in content:
                if p in sequences:
                    if sequences[p] in x_test:
                        duplicates.append(p)
                    if length_limit is None or len(sequences[p]) <= length_limit:
                        x_test.append(sequences[p])
                        y_test.append(y)

    return x_train, y_train, x_test, y_test
