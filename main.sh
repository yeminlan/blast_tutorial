#!/bin/bash

# 00.Create work directory
mkdir blast_tutorial
cd blast_tutorial


#1. Load BLAST module pre-installed on proteus
module list
module load ncbi-blast/gcc/64/2.2.29
blastn -h
blastn -help
# blastn/blastp/blastx/tblastn/tblastx
# makeblastdb -h


#2. Database to BLAST against
mkdir myblastdb
cd myblastdb

#2a. Build customized database from fasta file
#2a_1. Download nucleotides and proteins from an E.coli genome
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/Escherichia_coli_K_12_substr__W3110_uid161931/*.ffn
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/Escherichia_coli_K_12_substr__W3110_uid161931/*.faa

#2a_2. Build customized database from fasta files
makeblastdb -h
makeblastdb -in *.ffn -dbtype nucl -out Ecoli -parse_seqids
makeblastdb -in *.faa -dbtype prot -out Ecoli -parse_seqids
# "-parse_seqids" if the sequence headers were in NCBI style and therefore can be parsed. For example:
# ">gi|568336023|gb|CM000663.2| Homo sapiens chromosome 1, GRCh38 reference primary assembly."

#2b. Download existing database from ftp://ftp.ncbi.nlm.nih.gov/blast/db/
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz
tar zxvf swissprot.tar.gz
rm swissprot.tar.gz
cd ..


#3. What to BLAST?
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/Shigella_flexneri_5_8401_uid58583/*.ffn
# grep -n ">" *.ffn | head
head -n 91 *.ffn > Shigella.ffn

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/Shigella_flexneri_5_8401_uid58583/*.faa
#grep -n ">" *.faa | head
head -n 35 *.faa > Shigella.faa


#4. BLAST
blastn -db ./myblastdb/Ecoli -outfmt 6 -max_target_seqs 3 -query Shigella.ffn -out Shigella.fna.blast1
blastp -db ./myblastdb/Ecoli -outfmt 6 -max_target_seqs 3 -query Shigella.faa -out Shigella.faa.blast2
blastp -db ./myblastdb/swissprot -outfmt 6 -max_target_seqs 3 -num_threads 8 -query Shigella.faa -out Shigella.faa.blast3
blastn -remote -db nr -outfmt 6 -max_target_seqs 3 -query Shigella.ffn -out Shigella.fna.blast4
# -perc_identity 97
# -evalue 0.00001
# -best_hit_overhang 0.25 -best_hit_score_edge 0.25
# NOTE: use "-num_threads $NSLOTS" when submitting a proteus job script


#5. Parse BLAST output (using Shigella.faa.blast3 as an example)

# parse for hits matched to E.coli
grep "ECOLI" Shigella.faa.blast3

# parse for hits with >1000 Bit-score
awk '$12 > 1000 { print $0 }' Shigella.faa.blast3

# parse for best hits for each query
sort -u -k1,1 Shigella.faa.blast3 >  Shigella.faa.blast3.besthit

# count the hits by species
awk '{print $2}' Shigella.faa.blast3 | awk -F "_" '{print $2}' | sort | uniq -c

# count the hits by gene names
awk '{print $2}' Shigella.faa.blast3 | awk -F "|" '{print $5}' | awk -F "_" '{print $1}' | sort | uniq -c


#6. Tips for doing the same thing in python (optional)

# acquire sequences using BioPython
python biopython_get_seqs.py

# BLAST using BioPython
python biopython_blast.py

# parse BLAST output using BioPython
python biopython_blastparse_xml.py
python biopython_blastparse_tabular.py

