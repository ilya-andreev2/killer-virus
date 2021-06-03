## T27 - Phylogenetic Analysis of the DUP240 Gene Family

Gene trees were generated using *Fasttree 2.1* (double-precision version). Control files were prepared for *PAML 4.9j*. Alignment was generated using PRANK v.170427 using the following code:

```
prank -d=./RR_DUP240/REF_DUP240_subset.fsa -o=REF_DUP240_subset_curated.aln -codon -gaprate=0.000001 -gapext=0.001
```

Execute ```codeml <*.ctl>``` to perform analysis to test for presence of positive selection (models 7 vs 8). Use *subsetAlignment.jvp* (JalView) as a reference to back-calculate the site positions identified by BEB analysis.

Analysis excludes YHL044W, YCR007C, DUP-CLI3, and YAR023C as they are diverged from the genes included in the alignment, resulting in multiple unresolvable gaps.