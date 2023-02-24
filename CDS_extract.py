# CDS_extract.py
# Author: Shahd ElNaggar
# Date: 2.23.23
# This script finds all tetA and tetR genes in the plasmids based on representative BLAST hits, and aligns the sequences between all plasmids. 

from Bio import SearchIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import pandas as pd
import numpy as np
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO

# read files 
pth = "./tet_analysis/"
plas_file = pd.read_excel(pth + "plasmids_used.xlsx")
plas_list = plas_file[plas_file["drug"] == "tet"]["Plasmid"].values.tolist() # tet plasmid set 

if not os.path.exists(pth + "BLAST"): os.makedirs(pth + "BLAST")
if not os.path.exists(pth + "alignments"): os.makedirs(pth + "alignments")

files = [pth + 'plasmids/' + p + '/' + p + '.ffn' for p in plas_list] # multifasta DNA seq for each annotation
filesp = [pth + 'plasmids/' + p + '/' + p + '.faa' for p in plas_list] # protein sequence
locus_df = pd.DataFrame(plas_list)

## Direct BLAST Method to produce nucleotide sequences of genes (protein-protein blast)

# tetA, class C as query 
allrecs = []
tetA = []
for f,p,fn in zip(filesp,plas_list,files): #for each plasmid protein multi-fasta, blast against tetA class C
	outf = pth + "BLAST/curr.xml"
	blast_local = NcbiblastpCommandline(query = pth + "tetgenes/tetAC.txt", subject = f, evalue = 0.0005, outfmt = 5, out = outf)
	stdour,stderr = blast_local() 
	alns = [alignment.title.split(' ')[1] for record in NCBIXML.parse(open(outf,"r")) for alignment in record.alignments] 
	plas_records = SeqIO.to_dict(SeqIO.parse(fn,'fasta'))
	rec = [SeqRecord(Seq(plas_records[aln].seq), id = "{}_{}".format(p,aln), description = plas_records[aln].description) for aln in alns if aln in plas_records]
	allrecs += rec
	tetA += alns
SeqIO.write(allrecs,pth +"tetA_plasmids.fasta","fasta")
locus_df["tetA"] = tetA

# tetR as query
allrecs = []
tetR = []
for f,p,fn in zip(filesp,plas_list,files): #for each plasmid protein multi-fasta, blast against tetA class C
	outf = pth + "BLAST/curr.xml"
	blast_local = NcbiblastpCommandline(query = pth + "tetgenes/tetRD.txt", subject = f, evalue = 0.0005, outfmt = 5, out = outf) 
	stdour,stderr = blast_local() 
	alns = [alignment.title.split(' ')[1] for record in NCBIXML.parse(open(outf,"r")) for alignment in record.alignments] 
	plas_records = SeqIO.to_dict(SeqIO.parse(fn,'fasta'))
	rec = [SeqRecord(Seq(plas_records[aln].seq), id = "{}_{}".format(p,aln), description = plas_records[aln].description) for aln in alns if aln in plas_records]
	allrecs += rec
	tetR+= [alns]
SeqIO.write(allrecs,pth + "tetR_plasmids.fasta","fasta")
locus_df["tetR"] = tetR

locus_df.to_csv(pth + 'tet_locus.csv', index=False) # save locus tags 

## Alignment using ClustalW

# tetA
cline = ClustalwCommandline("clustalw2", infile= pth + "tetA_plasmids.fasta", output = "clustal", outfile = pth + "alignments/tetA_aln.aln") 
stdout, stderr = cline()
# tetR
cline = ClustalwCommandline("clustalw2", infile= pth + "tetR_plasmids.fasta", output = "clustal", outfile = pth + "alignments/tetR_aln.aln")
stdout, stderr = cline()




