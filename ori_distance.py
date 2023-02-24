# ori_distance.py
# Author: Shahd ElNaggar
# Date: 2.23.23
# This script finds the distance between the origin of replication/transfer and the tet gene of interest (tetA). 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

# oriC_tag
ori = "oriC"
oriT_tag = 0
if oriT_tag: ori = "oriT"

# read files
pth = "./tet_analysis/"
plas_file = pd.read_excel(pth + "plasmids_used.xlsx") 
plas_list = plas_file[plas_file["drug"] == "tet"]["Plasmid"].values.tolist() 
AC = plas_file[plas_file["drug"] == "tet"]["AC"].values.tolist()
gbks = [pth + 'plasmids/' + p + '/' + p + '.gbk' for p in plas_list] 
filesO = [pth + "{}_out/{}_{}.csv".format(ori,p,ori) for p in plas_list] 

lengths = []
for gb in gbks:
	for record in SeqIO.parse(gb,"genbank"):
		lengths.append(len(record.seq))

tetA_start = pd.read_csv(pth + "tetA_starts.csv") # tetA start positions

oriends = []
for o in filesO:
	oriCdf = pd.read_csv(o,names=["Name", "Locus", "Pct_ID", "3","4","5","oriC_start","oric_end","8","9","10","11"])
	end = oriCdf.iloc[[0],[7]]
	oriends.append(end.iloc[0,0])

tetA_start["ori_end"] = oriends
tetA_start["AC"] = AC 
tetA_start["Plasmid_Length"] = lengths
tetA_start["Distance_tetA_to_ori"] = tetA_start["Start"] + (tetA_start["Plasmid_Length"] - tetA_start["ori_end"])

tetA_start.to_csv(pth + "plasmid_distances_{}.csv".format(ori))

