# T27 - Phylogenetic Analysis of the DUP240 Gene Family

Gene trees were in the tree\ directory were generated using *Fasttree 2.1* (double-precision version; see [this blog post from the Darling lab](https://darlinglab.org/blog/2015/03/23/not-so-fast-fasttree.html)). Control files were prepared for *PAML 4.9j*. Alignment was generated using PRANK v.170427 using the following code:

```
prank -d=./REF_DUP240_subset.fsa -o=REF_DUP240_subset_curated.aln -codon -gaprate=0.000001 -gapext=0.001
```

GARD was run on [Datamonkey Adaptive Evolution Server](http://www.datamonkey.org/GARD/) and identified 7 potential recombination breakpoints, thus splitting everything into 8 segments (see log). Only files with prefix "7bp_" are used for analysis in the paper, with others included for completeness. 

Phobius predictions of transmembrane helices in DUPs is included

Execute ```codeml <*.ctl>``` to perform analysis to test for presence of positive selection (models 7 vs 8). Use *subsetAlignment.jvp* (JalView) as a reference to back-calculate the correct site positions identified by BEB analysis.

Analysis excludes YHL044W, YCR007C, DFP19, and YAR023C as they are too diverged from the genes included in the alignment.

Some of the files use older notation before DFP names were assigned, as follows:
* DUP YJM1 = DFP11
* DUP YJM2 = DFP12
* DUP YJM3 = DFP13
* DUP YJM4 = DFP14
* DUP YJM5 = DFP15
* DUP YJM6 = DFP16
* DUP CLI1 = DFP17
* DUP CLI2 = DFP18
* DUP CLI3 = DFP19
* DUP CLI4 = DFP20
* DUP CLI5 = DFP21
* DUP IFT1 = DFP22
* DUP IFT2 = DFP23
* DUP PWF1 = DFP24
* DUP YPS1 = DFP25
