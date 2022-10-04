#! /bin/bash

# prepare exonerate GFF output for Geneious
# need to pass experiment and paths to needed tar.gz files

usage () { echo "Usage : $0 -e <path to exonerate GFF> \
-p <path to 21295-exonerate_gff_to_alignment_gff3.pl perl file> \
-o <path to output GFF file> \

Specify the GFF file as output by exonerate.

Remove extraneous comment lines from GFF then modify CDS annotations to be GFF3 compliant.
"; }

echo $1

while getopts e:p:o: opt ; do
   case $opt in
	  e) ORIGINAL_GFF=$OPTARG ;;
	  p) PERL_CONVERSION=$OPTARG ;;
	  o) OUTPUT_GFF=$OPTARG ;;
	  *) usage; exit 1;;
   esac
done

if [ ! "$ORIGINAL_GFF" ] || [ ! "$PERL_CONVERSION" ] || [ ! "$OUTPUT_GFF" ] 
then
	usage
	exit
fi

 # generate new CDS annotations
 # save to temp file
 
new_cds_annotations_gff=$(mktemp /tmp/process-gff.XXXXXX)

perl "$PERL_CONVERSION" \
"$ORIGINAL_GFF" \
prot \
> "$new_cds_annotations_gff"

# remove header and footer
# then remove comment lines starting with #
# then remove original cds annotations
# save to temp file

cleaned_exonerate_gff=$(mktemp /tmp/process-gff.XXXXXX)

sed 's/-- completed exonerate analysis.*//' $ORIGINAL_GFF \
| sed 's/^#.*//' \
| sed 's/Command line:.*//' \
| sed 's/Hostname:.*//' \
 > "$cleaned_exonerate_gff"


# merge cleaned exonerate GFF and new CDS annotations
cat "$cleaned_exonerate_gff" "$new_cds_annotations_gff" \
| sed '/^\s*$/d' \
> "$OUTPUT_GFF"
