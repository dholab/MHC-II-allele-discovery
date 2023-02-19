#!/usr/bin/env python3

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, CompoundLocation

input_genbank = sys.argv[1]

ipd_seqs = SeqIO.parse(open(input_genbank),"genbank")

def removeSpecialCharacters(in_str, special_characters='*|: ', replace_character='_'):
    '''remove specified special characters from input str'''
    
    out_str = in_str.translate ({ord(c): replace_character for c in special_characters})
    
    return out_str

with open("ipd.fasta", "w") as ipd_output_handle:
    for record in ipd_seqs:
        ipd_record = SeqRecord(
        	Seq(str(record.seq).replace("X","")),
        	id=removeSpecialCharacters(record.name) + "_" + removeSpecialCharacters(record.id),
        	description=""
        )

        # add to reference
        # only add sequences longer than 100bp
        if len(record.seq) >= 100: SeqIO.write(ipd_record, ipd_output_handle, "fasta")