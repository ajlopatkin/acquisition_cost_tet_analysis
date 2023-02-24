#! /bin/bash

mkdir -p ./tet_analysis/oriT_out
mkdir -p ./tet_analysis/oriC_out

for filename in ./tet_analysis/plasmids/*/*.fna
do
	f=${filename##*/}
	f=${f%.*}
	blastn -query "${filename}" -subject ./tet_analysis/ori_db/orit.fna -out "./tet_analysis/oriT_out/${f}_oriT.csv" -outfmt 10
	blastn -query "${filename}" -subject ./tet_analysis/ori_db/oriC.fna -out "./tet_analysis/oriC_out/${f}_oriC.csv" -outfmt 10
done