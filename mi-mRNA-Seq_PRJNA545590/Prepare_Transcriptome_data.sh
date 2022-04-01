#!/bin/bash

#SBATCH --time=3:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=16   





mkdir Rat_data
cd Rat_data

#download Rat  transcriptome. All transcripts of Ensembl genes, excluding ncRNA.
wget http://ftp.ensembl.org/pub/release-105/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.mRatBN7.2.cdna.all.fa.gz


#download Rat Genome
wget http://ftp.ensembl.org/pub/release-105/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz

#unzip all
gunzip *.gz

#Make list of TranscriptIDs
grep ">" Rattus_norvegicus.mRatBN7.2.cdna.all.fa | awk '{print $1}' | sed 's/>//g' > ../TranscriptIDs


#create decoy list
grep ">" Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa | cut -d " " -f 1 | tr -d ">" > decoys.txt
        

#combine transcriptomes and decoy fasta files.   
cat  Rattus_norvegicus.mRatBN7.2.cdna.all.fa Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa > RatTranscrptome_decoy.fasta


#cleanup
rm Rattus_norvegicus.mRatBN7.2.cdna.all.fa Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa


#create salmon index
time salmon index -t RatTranscrptome_decoy.fasta -d decoys.txt -p 15 -i salmon_index

