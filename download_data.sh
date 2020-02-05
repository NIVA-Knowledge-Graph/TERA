#!/bin/bash

mkdir tmp

wget ftp://newftp.epa.gov/ecotox/ecotox_ascii_12_12_2019.exe
unzip ecotox_ascii_12_12_2019.exe
rm ecotox_ascii_12_12_2019.exe
mv ecotox_ascii_12_12_2019 ecotox_data

#enter input encoding here
FROM_ENCODING="Windows-1252"
#output encoding(UTF-8)
TO_ENCODING="UTF-8"
#convert
CONVERT=" iconv  -f   $FROM_ENCODING  -t   $TO_ENCODING"
#loop to convert multiple files 
for  file  in ./ecotox_data/*.txt; do
     $CONVERT   "$file"   -o  "${file%.txt}.txt"
done
for  file  in ./ecotox_data/validation/*.txt; do
     $CONVERT   "$file"   -o  "${file%.txt}.txt"
done

wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip --directory-prefix=tmp/
unzip tmp/new_taxdump.zip -d ./taxdump
rm -r tmp

mkdir pubchem
wget ftp://ftp.ncbi.nlm.nih.gov/pubchem/RDF/compound/general/pc_compound2parent.ttl.gz --directory-prefix=pubchem/
wget ftp://ftp.ncbi.nlm.nih.gov/pubchem/RDF/compound/general/pc_compound_type.ttl.gz --directory-prefix=pubchem/
gzip -d pubchem/*

mkdir chemble
wget ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBL-RDF/latest/* --directory-prefix=chemble/
gzip -d chemble/*

mkdir mesh
wget ftp://nlmpubs.nlm.nih.gov/online/mesh/rdf/2019/mesh2019.nt.gz --directory-prefix=mesh
gzip -d mesh/*

mkdir eol
wget https://editors.eol.org/other_files/SDR/traits-20190822/traits_20190822.zip --directory-prefix=eol/
unzip eol/*
