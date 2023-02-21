Hedgehog protein family in Serpentes

The hedgehog protein family is important for cell-cell signalling to direct cell differentiation during embryonic development in multicellular organisms. Snakes have three Hedgehog paralogs- Sonic, Desert, and Indian. A database was compiled of 8 Serpentes species each with all three paralogs. Sequences were found on NCBI GenBank.

Steps for data:
-Searched for python bivitattus hedgehog protein on NCBI GenBank
-Performed BLASTP on sequence against Serpentes
-Compiled nucleotide and amino acid information for blast results, excluding sequences labeled as "hypothetical" and species without all three paralogs

MSA with ClustalW
ClustalW2 was chosen to align because it is currently the only MSA option which I've been able to run. ClustalW is limited by the fact it is a progressive alignment which builds from the leaves. Nodes cannot be rearranged following alignment because it assumes the most probable alignment of the first two sequences is true in relation to all other sequences.
-copied aa sequences to txt file then converted it to fasta file
-stored file in bot563/ClustalW2
-in ClustalW2 program, data file loaded in using option 1 "Sequence Input From Drive"
-chose option 2 "Multiple Alignments"
-chose option 1 "Do complete multiple alignment now"
