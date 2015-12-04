Sequence alignment statistics were calculated using [AMAS](https://peerj.com/preprints/1355/). The version used is v0.95, Git commit [84a679a](https://github.com/marekborowiec/AMAS/commit/84a679ac71bc64ef94cb9d606dd535ae82226e25).
The statistics based on trees or on trees and sequences (regression) were calculated in R using `ape` (Paradis 2012) and `seqinr` (Charif et al. 2007) libraries. The code used is in `tree_props.R`. It is a slightly modified version of the code for the metazoan phylogeny paper (Borowiec 2015; [GitHub](https://github.com/marekborowiec/metazoan_phylogenomics/blob/master/gene_stats.R)) with “clock-likeness” measure added.
Average sequence heterogeneity was calculated using p4 (Foster 2004). This is the `aln_hetero.py` script.
A brief explanation of what is being calculated and how in the file `good_gene_stats.csv`; For each alignment/corresponding tree I calculated:
* alignment length
* number of taxa
* total matrix cells
* count of undetermined characters (X, N, O, -, ?)
* percent missing (using the above)
* number of variable sites (undetermined chars are excluded when determining this)
* proportion of variable sites
* number of parsimony informative sites (excluding undetermined)
* proportion of parsimony informative sites
* AT content
* GC content
* counts of all nucleotides, gap -, and missing ?
* average matrix heterogeneity (this is mean of Euclidean distance matrix of compositions)
* average bootstrap
* average branch length (including internal branches; possible rate proxy, lower number means slower-evolving)
* “clocklikeness” score (this is a measure how close to ultrametric a tree is; the algorithm finds a root that minimizes coefficient of variation in root to tip distances and returns that value; lower value is more clock-like: ultrametric tree has a score of 0)
* average uncorrected p-distance
* regression slope of identity distances plotted against branch lengths (the higher the value the closer the alignment is fitting to linear regression, which means lower saturation potential)
* R-squared of regression (as above, higher means better fit to linear regression and less saturation potential)

To cite:

Borowiec, M. L., Lee, E. K., Chiu, J. C., & Plachetzki, D. C. 2015. Extracting phylogenetic signal and accounting for bias in whole-genome data sets supports the Ctenophora as sister to remaining Metazoa. BMC Genomics; 16(1):987.

Charif D, Lobry JR. 2007. SeqinR 1.0-2: a contributed package to the R project for statistical computing devoted to biological sequences retrieval and analysis in Structural approaches to sequence evolution: Molecules, networks, populations (U. Bastolla, M. Porto, H.E. Roman and M. Vendruscolo Eds.) Biological and Medical Physics, Biomedical Engineering; pp 207–232.

Foster, P. G. 2004. Modeling compositional heterogeneity. Systematic Biology; 53(3):485–495.
Paradis E, Claude J, Strimmer K. 2004. APE: Analyses of phylogenetics and evolution in R language. Bioinformatics; 20:289–90.
