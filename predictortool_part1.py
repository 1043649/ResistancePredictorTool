# This code imports a fasta file of nucleotide sequences and translates them, 
# then calculates the most common sequence length 
# and uses sequences of this length to define a set of motifs



pip install biopython 
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import statistics



#### GENERATING AMINO ACID SEQUENCES ####

# import the fasta file of nucleotide sequences containing all the unique alleles of gyrA, as the variable 'ntsequences':
ntsequences = SeqIO.parse("ResistancePredictorTool/ExampleData/gyrAuniquealleles.fas", "fasta")

# write a new file, AAseqs.fasta, which will contain the translated amino acid sequneces:
with open("AAseqs.fas", "w") as aa_fas:
    for dna_record in ntsequences:
        dna_seqs = [dna_record.seq]
        # generate all translation frames, use bacterial codon table (NCBI table 11) and stop translating if a stop codon is reached:
        aa_seqs = (dnaseq[i:].translate(table=11, to_stop=True) for i in range(3) for dnaseq in dna_seqs)
        # select the sequence from the correct translation frame (the longest sequence):
        max_aa = max(aa_seqs, key=len)
        # write a new sequence record in the file 'AAseqs.fasta':
        aa_record = SeqRecord(max_aa, id=dna_record.id, description="translated sequence")
        SeqIO.write(aa_record, aa_fas, "fasta")
        
# import the file of translated sequences as the variable 'aasequences':
aasequences = SeqIO.parse("/Users/1043649/predictortoolcode/AAseqs.fas", "fasta")

# extract the lengths of the sequences:
seq_lengths = []
for sequence in list(aasequences):
    seq_lengths.append(len(sequence))

# find the most common sequence length:
modeseqlen = statistics.mode(seq_lengths)
print("Most common sequence length: ", modeseqlen)



#### USER INPUTS ####
# ask the user for their input (for gyrA inputs are: position = 91, WT = S, resistance mutations = F Y I, other mutations/ unknown = T):

position_index = (int(input("What is the position of the mutation?")) -1) # position of the mutation (store as the index of that position)
wt_aa = input("What is the amino acid in the wild type?") # wild type amino acid  

mut_aa_input = input("What are the mutated amino acids which predict the phenotype of interest?") # mutated amino acid (you can input multiple letters here by adding a space between them)
mut_aa = mut_aa_input.split() # creates a list if you have inputted multiple letters

other_aa_input = input("Are there any mutations with an unknown phenotype?") # other amino acid(s) with an unknown phenotype
other_aa = other_aa_input.split()

out_wt = input("What is the phenotype when the WT amino acid is present?") #ie. sensitive
out_mut = input("What is the phenotype when the mutated amino acid is present?") #ie. resisant
out_other = "unknown" 

# print the inputs to check they are correct:
print("Index of position: " + str(position_index))
print("WT amino acid: " + wt_aa)
print("Mutated amino acid: ", mut_aa)
print("Amino acid(s) with unknown phenotype: ", other_aa)
print("Outputs: " + out_wt + ", " + out_mut + ", " + out_other)



#### DEFINING THE MOTIFS ####

# import the file of the unique amino acid sequences created earlier
aasequences = SeqIO.parse("/Users/1043649/predictortoolcode/AAseqs.fas", "fasta")

# extract only the sequences which have the most common length (to account for any sequencing issues which have led to slightly longer or shorter sequences):
motif_sequences = []
for sequence in list(aasequences):
    if len(sequence) == modeseqlen:
        motif_sequences.append(sequence.seq)
print("Number of sequences with length " + str(modeseqlen) + ": " + str(len(motif_sequences)))

# create a new list that is populated by the sequences which have either the WT or the mutated amino acids at the specified position:
motif_sequences_final = []
for sequence in motif_sequences:
    if sequence[position_index] == wt_aa:
        motif_sequences_final.append(sequence)
    for item in mut_aa:
        if sequence[position_index] == item:
            motif_sequences_final.append(sequence)
    for item in other_aa:
        if sequence[position_index] == item:
            motif_sequences_final.append(sequence)
print("Number of sequences with known amino acids at position " + str((position_index + 1)) + ": " + str(len(motif_sequences_final)))

# extract the 10 amino acids from either side of the position to define the 21-amino acid motifs:
motifs = []
for seq in motif_sequences_final:
    motifs.append(seq[(position_index - 10):(position_index + 11)])
motifs_unique = []
for seq in motifs:
    if seq not in motifs_unique:
        motifs_unique.append(seq)
print("Unique motifs: ", motifs_unique)
print("Number of unique motifs: ", len(motifs_unique))

# identify which of the motifs predict which phenotype (ie. which are sensitive vs resistant):
for sequence in motifs_unique:
    if sequence[10] == wt_aa:
        print(sequence + " = " + out_wt)
    elif sequence[10] in mut_aa:
        print(sequence + " = " + out_mut)
    else:
        print(sequence + " = unknown")


        
  
#### CREATE A FASTA FILE CONTAINING THE MOTIFS ####  

# calculate how many of the motifs are sensitive/resistant/unknown 
num_sensitive = 0
num_resistant = 0
num_unknown = 0
for sequence in motifs_unique:
    if sequence[10] == wt_aa:
        num_sensitive += 1
    elif sequence[10] in mut_aa:
        num_resistant += 1
    else:
        num_unknown += 1
print(num_sensitive, num_resistant, num_unknown)

# make a list of the motif sequence names for the fasta file:
motifnames_resistant = []
for i in range(1, num_resistant+1):
    motifnames_resistant.append(">" + out_mut + str(i))
motifnames_sensitive = []
for i in range(1, num_sensitive+1):
    motifnames_sensitive.append(">" + out_wt + str(i))
motifnames_unknown = []
for i in range (1, num_unknown+1):
    motifnames_unknown.append(">" + out_other + str(i))
print(motifnames_resistant)
print(motifnames_sensitive)
print(motifnames_unknown)

resistant_motifs = []
for sequence in motifs_unique:
    if sequence[10] in mut_aa:
        resistant_motifs.append(sequence)
        
sensitive_motifs = []
for sequence in motifs_unique:
    if sequence[10] == wt_aa:
        sensitive_motifs.append(sequence)

unknown_motifs = []
for sequence in motifs_unique:
    if sequence[10] in other_aa:
        unknown_motifs.append(sequence)

# create a dictionary containing all the motif sequences and sequence names:        
resistant_dict = dict(zip(motifnames_resistant, resistant_motifs))
sensitive_dict = dict(zip(motifnames_sensitive, sensitive_motifs))
unknown_dict = dict(zip(motifnames_unknown, unknown_motifs))
motifs_dict = {**resistant_dict, **sensitive_dict, **unknown_dict}

# create a fasta file using the motifs using this dictionary:
with open("motifs.fas", "w") as f:
    for key, value in motifs_dict.items():
        f.write("%s\n%s\n" % (key, value))
