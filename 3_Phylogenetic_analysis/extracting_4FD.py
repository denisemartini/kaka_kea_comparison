import sys
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.Data.CodonTable import unambiguous_dna_by_id
from Bio.Data import CodonTable
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

# 23.01.19 Denise Martini
# This script takes an alignment file in clustal format (output from tcoffee)
# and extracts 4-fold degenerate sites.
# The argument to feed to the script is the partial name of the file
# and the file needs to end in .best.fas
# e.g. for a file named TF_1047.best.fas your run the script like:
# python extracting_4FD.py TF_1047

# reading in the argument (change this if your file ends in a different suffix)
infile = sys.argv[1] + ".best.fas"
print(infile)

# importing the alignment
alignment = AlignIO.read(infile, "fasta")
print("Alignment length %i" % alignment.get_alignment_length())

# defining the functions to find degenerate codons
# (from https://biopython.org/wiki/Degenerated_Codons)


def altcodons(codon, table):
    """List codons that code for the same aminonacid / are also stop.

    @param codon
    @table code table id
    @return list of codons

    """
    tab = unambiguous_dna_by_id[table]

    if codon in tab.stop_codons:
        return tab.stop_codons

    try:
        aa = tab.forward_table[codon]
    except:
        return []

    return [k for (k, v) in tab.forward_table.items()
            if v == aa and k[0] == codon[0] and k[1] == codon[1]]


def degeneration(codon, table):
    """Determine how many codons code for the same amino acid / are also stop

    @param codon the codon
    @param table code table id
    @param the number of codons also coding for the amino acid codon codes for

    """
    return len(altcodons(codon, table))


def is_x_degenerated(x, codon, table):
    """Determine if codon is x-fold degenerated.

    @param codon the codon
    @param table code table id
    @param true if x <= the degeneration of the codon

    """
    return (x <= len(altcodons(codon, table)))


def degenerated_subseq(seq, x, table):
    """Get a subsequence consisting of the x-fold degenerated codons only."""
    data = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        if is_x_degenerated(x, codon, table):
            data.extend([i])
    return data


# identifying 4-fold degenerate codons in the every sequence of the alignment
deg = []
for record in alignment:
    deg.append(degenerated_subseq(str(record.seq), 4, 1))
# intersecting the results from all sequences
all_deg = sorted(set(deg[0]).intersection(*deg[:len(deg)]))

# selecting the third position from each codon identified
final_deg = []
for n in all_deg:
    final_deg.append(n + 2)

# extracting the identified 4FD sites from each sequence in the alignment
seqs = []
for record in alignment:
    edited = []
    for n in final_deg:
        edited.append(str(record.seq[n]))
    seqs.append(''.join(edited))

# merging them in a new alignment, retrieving IDs from original file
final_alg = MultipleSeqAlignment([
    SeqRecord(Seq(str(seqs[0]), generic_dna), id=str(alignment[0].id)),
    SeqRecord(Seq(str(seqs[1]), generic_dna), id=str(alignment[1].id)),
    SeqRecord(Seq(str(seqs[2]), generic_dna), id=str(alignment[2].id)),
    SeqRecord(Seq(str(seqs[3]), generic_dna), id=str(alignment[3].id)),
    SeqRecord(Seq(str(seqs[4]), generic_dna), id=str(alignment[4].id)),
    SeqRecord(Seq(str(seqs[5]), generic_dna), id=str(alignment[5].id)),
         ])

# defining the name of the output file, change here for a different suffix
filename = sys.argv[1] + "_nt_4FD.aln"
print(filename)

# outputting the edited alignment in fasta format
AlignIO.write([final_alg], filename, "fasta")
