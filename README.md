# acquisition_cost_tet_analysis
Analysis scripts for tetracycline operon 

readMe.txt
author: Shahd ElNaggar

1. plasmids_used.xlsx contains the plasmid metadata.
2. CDS_extract.py BLASTp plasmid proteins against representative tetA and tetR sequences; saves tetA/R nucleotide sequences in multi-fasta; aligns sequences using ClustalW; generates list of tetA/R locus tags in genbank files. 
2. tet_operon.py uses locus tags generated by step (1) to extract nucleotide sequence of tetracycline operon from genbank files; it also creates a file containing tetA start positions.
3. ori_BLAST.sh to BLAST plasmids against oriC/T databases.
4. oriC_tet.py calculates distance between origin of transfer/replication and tetA gene to create a csv file with distances for each gene. 

dependencies: BioPython, ClustalW, BLAST
