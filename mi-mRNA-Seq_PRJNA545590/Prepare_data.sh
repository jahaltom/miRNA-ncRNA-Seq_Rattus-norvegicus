#!/bin/bash

#SBATCH --time=3:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=16   


#Download Rat gtf
wget http://ftp.ensembl.org/pub/release-105/gtf/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.105.gtf.gz

mkdir Rat_data
cd Rat_data

#download Rat ncrna transcriptome 
wget http://ftp.ensembl.org/pub/release-105/fasta/rattus_norvegicus/ncrna/Rattus_norvegicus.mRatBN7.2.ncrna.fa.gz


#download Rat Genome
wget http://ftp.ensembl.org/pub/release-105/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz

#unzip all
gunzip *.gz



#create decoy list
grep ">" Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa | cut -d " " -f 1 | tr -d ">" > decoys.txt
        

#combine transcriptomes and decoy fasta files.   
cat  Rattus_norvegicus.mRatBN7.2.ncrna.fa Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa > RatTranscrptome_decoy.fasta


#cleanup
rm Rattus_norvegicus.mRatBN7.2.ncrna.fa Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa


#create salmon index
time salmon index -t RatTranscrptome_decoy.fasta -d decoys.txt -p 15 -i salmon_index

