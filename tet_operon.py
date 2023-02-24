# tet_operon.py
# Author: Shahd ElNaggar
# Date: 2.23.23
# Extracting entire nucleotide sequence from tetR to tetA in plasmids

from Bio import SearchIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
import pandas as pd
import numpy as np

pth = "./tet_analysis/"
plas_file = pd.read_excel(pth + "plasmids_used.xlsx")
plas_list = plas_file[plas_file["drug"] == "tet"]["Plasmid"].values.tolist() #list of plasmids with tet data 
gbks = [pth + 'plasmids/' + p + '/' + p + '.gbk' for p in plas_list] 

tet_loci = pd.read_csv(pth + "tet_locus.csv")

# formatting
tetR = tet_loci["tetR"].str.strip("',[]")
tet_loci["Loci"] = tet_loci["tetA"] + " " + tetR.str.replace("'","")
tet_loci["Loci"] = tet_loci["Loci"].str.replace(",","")
df = tet_loci[["0","Loci"]]

df2 = pd.DataFrame(index = plas_list) 

i = 0 
allrecs = []
tetA_starts = []
for p,f in zip(plas_list,gbks): 
	loci = df["Loci"][i]
	i += 1
	loci = loci.split(" ")
	loci.sort()
	curr_record = SeqIO.parse(f,"genbank")
	for record in curr_record: 
		for feature in record.features:
			if feature.type == "CDS":
				tag = feature.qualifiers.get('locus_tag')[0]
				if tag in loci: 
					if "gene" in feature.qualifiers: # tetA start position
						if feature.qualifiers['gene'][0] == "tetA":
							tetA_start = feature.location._start.position
							tetA_starts.append(tetA_start)
					if loci.index(tag) == 0: 
						start = (feature.location._start.position)
					elif loci.index(tag) == (len(loci)-1):
						end = (feature.location._end.position)
		myseq = SeqRecord(Seq(record.seq[start:end]),id = str(p), description = "Tet Operon Estimate")
		allrecs.append(myseq)
SeqIO.write(allrecs,"tet_analysis/tetoperonestimates.fasta","fasta")
df2["Start"] = tetA_starts
df2.to_csv(pth + "tetA_starts.csv")

## CLUSTALW for operon 
cline = ClustalwCommandline("clustalw2", infile= pth + "tetoperonestimates.fasta", output = "clustal", outfile = pth + "alignments/tetoperonestimate.aln") 
stdout, stderr = cline()


