# if header (starts with @) then simply print line
# if not header, replace 'N' in CIGAR ($6) with 'D'
BEGIN{OFS="\t"}
{
if ($1 ~ /^@/) 
{ 
  print $0; 
}
else 
{
  gsub("N", "D", $6); print $0
}
}