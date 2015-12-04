#! /usr/bin/env python

from __future__ import division
import glob, p4

# stop ckecking for identical sequences within alignments
p4.var.doCheckForDuplicateSequences = False

# get all alignments
alignments = glob.glob('./numbered_alignments/*fas')

# calculate mean
def mean(a):
    return sum(a) / len(a)

# calculate mean compositional distance
def composition(file_name):
    alignment = p4.func.readAndPop(file_name)
    comp_matrix = alignment.compositionEuclideanDistanceMatrix()
    mean_comp = mean(map(mean, zip(*comp_matrix.matrix)))
    return mean_comp

# create a dict with alignment names as keys and mean comp dist as values
aln_dict = {aln:composition(aln) for aln in alignments}

# write to a file
out_file_name = "heterogeneity.txt"

out_file = open(out_file_name, "w")

for file_name, mean_comp in aln_dict.items():
    out_file.write('{},{}'.format(file_name, mean_comp))

out_file.close()
