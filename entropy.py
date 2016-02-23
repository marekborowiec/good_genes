#! /usr/bin/env python3

import math
from collections import Counter
import numpy as np
from amas import AMAS
from sys import argv

np.set_printoptions(threshold=np.nan)

in_f = argv[1]
frmt = argv[2]

meta_aln = AMAS.MetaAlignment(in_files=[in_f], data_type="dna",in_format=frmt, cores=1)

chars = 'ACGT'

aln_dicts = meta_aln.get_parsed_alignments()
aln = aln_dicts[0]
aln_len = len(list(aln.values())[0])

matrix = [list(sequence) for sequence in aln.values()]

def entropy_calc(p):
    p = p[p != 0]
    return np.dot(-p,np.log2(p)) # the function returns the entropy resul

def get_column(i):
    return [row[i] for row in matrix]

def get_no_missing_ambig(matrix):
    no_miss_amb = []
    for column in range(aln_len):
        site = get_column(column)
        new_site = [char for char in site if char in chars]
        no_miss_amb.append(new_site)
    return no_miss_amb

new_matrix = get_no_missing_ambig(matrix)

all_char_counts = []
for column in new_matrix:
    column_counts = []
    for char in chars:
        count = column.count(char)
        column_counts.append(count)
    all_char_counts.append(column_counts)

arr = np.array(all_char_counts)

#print(arr)
props = arr/arr.sum(axis=1, keepdims=True)

#print(props)
site_entropy = np.array([[entropy_calc(site)] for site in props])
site_entropy[np.isnan(site_entropy)] = 0

pos = 'ABC'
#for index, site in enumerate(site_entropy):
 #   print('site{}{},{}'.format((index + 1), pos[int(index % 3)], float(site)))

print('{}\t{}'.format(in_f, (np.sum(site_entropy) / aln_len)))
#print(aln_len)
