#!/usr/bin/env python3
from Bio import SeqIO
import sys

def removeSpecialCharacters(in_str, special_characters='*|: ', replace_character='_'):
    '''remove specified special characters from input str'''
    
    out_str = in_str.translate ({ord(c): replace_character for c in special_characters})
    
    return out_str

# command line parameter of FASTA sequences to process
sequences_to_rename = sys.argv[1]

# Open the input and output files
with open(sequences_to_rename, 'r') as input_file, open('all_fasta_renamed.fasta', 'w') as output_file:
    
    # Parse the input file and write the updated records to the output file
    records = SeqIO.parse(input_file, 'fasta')

    for record in records:
        record.name = removeSpecialCharacters(record.name)
        record.id = removeSpecialCharacters(record.id)
        SeqIO.write(record, output_file, 'fasta')