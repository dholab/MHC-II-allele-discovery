#!/usr/bin/env python3

import sys
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

cdna_matches = sys.argv[1]
putative_animal = sys.argv[2]
cdna_ref = sys.argv[3]

# create gDNA sequence object
for gdna_record in SeqIO.parse(cdna_matches, "fasta"):
    # print status message
    print("Finding cDNA matches to " + str(gdna_record.name))
    
    # map gDNA to cDNA
    # returning FASTA file of cDNA within 50nt
    # write gDNA sequence to file
    with open(putative_animal + "_gdna_single_temp.fasta", "w") as output_handle:
        SeqIO.write(gdna_record, output_handle, "fasta")
        
    # create cDNA sequence object from nearby matches
    for cdna_record in SeqIO.parse(cdna_ref, "fasta"):

        # write a tab delimited line containing the gdna_record.id, cdna_record.id, cdna_record length, and match length
        # use command substitution to get the match length
        # run MUSCLE on sequence objects, but do it in a stream to avoid a lot of file I/O
        # pipe output to CLUSTALW format
        # then count the number of '*' characters that denote matches between the two sequences
        # this works for class I
        # for class II, the cDNA can be longer than the gDNA so this doesn't work
        # if the count of natching characters equals the number of cDNA characters, write to file
        # add maxiters = 2 to accelerate processing per https://www.drive5.com/muscle/manual/compromise.html
        
        merged_output = putative_animal + "_merged.aln"
        
        command = 'echo "' + gdna_record.name  + '\t' + cdna_record.name + '\t' + str(len(gdna_record.seq)) + '\t' + str(len(cdna_record.seq)) + '\t' + ' \
            $(echo ">' + gdna_record.name + '\n' + str(gdna_record.seq) + '\n>' + cdna_record.name + '\n' + str(cdna_record.seq) + '" \
            | muscle -maxiters 2 -quiet -clwstrict \
        | grep "^ " | grep "\*" -o | wc -l )"' + '>>' + merged_output
        
        subprocess.call(command, shell=True)
