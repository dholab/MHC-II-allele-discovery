#!/usr/bin/env python3

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, CompoundLocation

input_fasta = sys.argv[1]

ipd_seqs = SeqIO.parse(open(input_fasta),"genbank")

def removeSpecialCharacters(in_str, special_characters='*|: ', replace_character='_'):
    '''remove specified special characters from input str'''
    
    out_str = in_str.translate ({ord(c): replace_character for c in special_characters})
    
    return out_str

with open("gdna_reference.fasta", "w") as gdna_output_handle:
  with open("cdna_reference.fasta", "w") as cdna_output_handle:
    for record in ipd_seqs:
      
      # Extract exon locations from the file
      exon_locations = []
      for feature in record.features:
          if feature.type == "exon":
              exon_location = FeatureLocation(start=feature.location.start, end=feature.location.end)
              exon_locations.append(exon_location)

      # if there are no exon_locations, skip sequence
      if len(exon_locations) == 0:
          print(record.name)
          continue

      # if there is only a single exon location
      # it is an exon 2 sequence
      # this should be added to the gDNA reference FASTA
      if len(exon_locations) == 1:
        print(record.name)
        # create gdna record
        gdna_record = SeqRecord(
        	Seq(str(record.seq).replace("X","")),
        	id=removeSpecialCharacters(record.name) + "_" + removeSpecialCharacters(record.id),
        	description=""
        )

        # add to reference
        # only add sequences longer than 100bp
        if len(record.seq) >= 100: SeqIO.write(gdna_record, gdna_output_handle, "fasta")

          
      # make compound location
      # test if length of sequence equals length of compound location
      # if so, exons are joined and it is a cDNA sequence
      # if not, there is sequence between exons and it is a gDNA sequence
      if len(exon_locations) > 1:
        compound_location = CompoundLocation(exon_locations)
        
        # if length of entire sequence equals the length of the compound_location
        # is indicates concatenated exons -> cDNA sequence
        # otherwise, it is a gDNA sequence
        if len(compound_location) == len(record):
          # cDNA
          # create cdna record
          cdna_record = SeqRecord(
            Seq(str(record.seq).replace("X","")),
            id=removeSpecialCharacters(record.name) + "_" + removeSpecialCharacters(record.id),
            description=""
          )
          if len(record.seq) >= 100: SeqIO.write(cdna_record, cdna_output_handle, "fasta")
        
        else:
          # gDNA
          # create gdna record
          gdna_record = SeqRecord(
            Seq(str(record.seq).replace("X","")),
            id=removeSpecialCharacters(record.name) + "_" + removeSpecialCharacters(record.id),
            description=""
          )      
          if len(record.seq) >= 100: SeqIO.write(gdna_record, gdna_output_handle, "fasta")
