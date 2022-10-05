#!/usr/bin/env python3

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

input_fasta = sys.argv[1]
animal = sys.argv[2]

ipd_seqs = SeqIO.parse(open(input_fasta),"genbank")

def removeSpecialCharacters(in_str, special_characters='*|: ', replace_character='_'):
    '''remove specified special characters from input str'''
    
    out_str = in_str.translate ({ord(c): replace_character for c in special_characters})
    
    return out_str

with open(str(animal) + "_gdna_reference.fasta", "w") as gdna_output_handle:
  with open(str(animal) + "_cdna_reference.fasta", "w") as cdna_output_handle:
    for record in ipd_seqs:
      # count number of intron features
      intron_ct = 0
      if record.features:
        for feature in record.features:
              # if intron features, write to gDNA, otherwise write to cDNA
              if (feature.type == "intron"):
                intron_ct += 1

      # if intron_ct > 0 then gDNA
      if intron_ct == 0: 
        # create gdna record
        gdna_record = SeqRecord(
        	Seq(str(record.seq).replace("X","")),
        	id=removeSpecialCharacters(record.name),
        	description=""
        )

        # only write records at least 100nt
        if len(record.seq) >= 100: SeqIO.write(gdna_record, cdna_output_handle, "fasta")
      else:
        # create cdna record
        cdna_record = SeqRecord(
        	Seq(str(record.seq).replace("X","")),
        	id=removeSpecialCharacters(record.name),
        	description=""
        )
        
        # only write records at least 100nt
        if len(record.seq) >= 100: SeqIO.write(cdna_record, gdna_output_handle, "fasta")
