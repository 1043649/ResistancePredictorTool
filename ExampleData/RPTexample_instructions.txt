To run the resistance predictor tool (using the Neisseria gonorrhoeae GyrA_91 mutation to predict ciprofloxacin resistance as an example) see the following instructions. Additional downloads/applications required are: Python and a Python package installer such as Anaconda, a Python interpreter application (eg. Jupyter Lab), and NCBI BLAST+. 


1. In a Python interpreter, run the script 'predictortool_part1.py'. Use the file 'gyrAuniquealleles.fas' as input: this contains all the unique gyrA nucleotide sequences recorded in PubMLST as of 12th May 2023, which can be downloaded from https://pubmlst.org/bigsdb?db=pubmlst_neisseria_seqdef&page=alleleQuery&locus=NEIS1320. The user inputs for GyrA_91 are: position = 91, WT = S, resistance mutations = F Y I, other mutations/ unknown = T).

2. The database of motifs created when running part 1 of the tool is available as the file 'gyrA91motifs.fas'.

3. In Terminal (Mac) or Command-line (Windows), run the first command in 'predictortool_part2.txt' to generate the BLAST database using these motifs. Then run the second command to search a query sequence against this BLAST database. An example query sequence is available as the file 'WHOYgyrA.fas': this is the gyrA sequence from the Neisseria gonorrhoeae WHO reference strain Y (Unemo, et al., 2016). 

4. Identify the predicted phenotype from the resulting output file (this file is available as 'WHOYblastx.txt'). The phenotype of the motif which is the top BLAST hit predicts the phenotype of the query sequence. In this example, WHO-Y is predicted to be resistant which corresponds to the laboratory-reported ciprofloxacin SIR of this strain.  
