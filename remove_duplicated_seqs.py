#!/usr/bin/python
# ***************************************************************
# Name:      remove_duplicated_seqs.py
# Purpose:   This script removes duplicated sequences in fasta files
#            according to sequence (rather than header), and outputs
#            the duplicated sequences with the matching sequence IDs
#            in deduplicated sequence file.
# Ref: https://biopython.org/wiki/Sequence_Cleaner
# Dependencies:
#            Biopython is installed
# Version:   0.1
# Authors:   Zihan DAI (z.dai.1@research.gla.ac.uk)
# Created:   2019-02-14
# ***************************************************************

import sys
import argparse
import Bio

import warnings
from Bio import BiopythonExperimentalWarning
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)
    from Bio import SeqIO

def parse_args():
	parser = argparse.ArgumentParser(description = '%(prog)s removes duplicated sequences in fasta files according to sequence (rather than header), and outputs the duplicated sequences with the matching sequence IDs in deduplicated sequence file.')
	parser.add_argument('-i', '--input-file', dest='fasta_file', help = 'input fasta file', required = True)
	return parser.parse_args()

def add_duplicated_seq_ID(seq, id, dupe_sequences, uniq_sequences):
    # Check if the sequence is already in dupe_sequences
    if seq not in dupe_sequences:
        dupe_sequences[seq] = uniq_sequences[seq] + "," + id
    else:
        dupe_sequences[seq] += "," + id
    return dupe_sequences

def remove_dupe(fasta_file, dupe_sequences, uniq_sequences):
    # Read fasta input
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        # Take the current sequence
        seq = str(seq_record.seq).upper()
        # Check if the sequence is already in the dictionary
        if seq not in uniq_sequences:
            uniq_sequences[seq] = seq_record.id
        # If the sequence is already in the dictionary, the duplicated sequence
        # would be added into dupe_sequences as key, and
        else:
            add_duplicated_seq_ID(seq, seq_record.id, dupe_sequences, uniq_sequences)
    return uniq_sequences, dupe_sequences

def main():
    args = parse_args()
    fasta_file = args.fasta_file

    # Create dictionary to add the sequences
    uniq_sequences={}
    # Create dictionary to add the removed sequence IDs
    dupe_sequences={}

    remove_dupe(fasta_file, dupe_sequences, uniq_sequences)

    # Write the outputs
    with open("unique_" + fasta_file, "w+") as output_file_1:
        # Just read the hash table and write on the file as a fasta format
        for seq in uniq_sequences.keys():
            output_file_1.write(">" + uniq_sequences[seq] + "\n" + seq + "\n")

    with open("duplicated_" + fasta_file, "w+") as output_file_2:
        for seq in dupe_sequences.keys():
            output_file_2.write(">" + dupe_sequences[seq] + "\n" + seq + "\n")

if __name__ == '__main__':
	main()
