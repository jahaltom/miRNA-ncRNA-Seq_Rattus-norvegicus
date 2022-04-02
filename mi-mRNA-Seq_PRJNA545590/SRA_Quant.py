import yaml
import sys
import os
from pyrpipe import sra,quant,qc
import pandas as pd


#Read in config file
configfile: "config.yaml"
DIR = config['DIR']
THREADS=config['THREADS']
Tr=config['Transcriptome']


#creates a trim_galore object.
trim_galore=qc.Trimgalore(threads=4)

#Read in run_accession ids
with open ("RAids.txt") as f:
        ra_ids=f.read().splitlines()

#create salmon object. If index is not present it will be created
salmon=quant.Salmon(index=Tr.split("/")[0]+"/salmon_index",transcriptome=Tr,threads=THREADS)

rule all:
        input:
                expand("{wd}/results_TPM_tx.tsv",wd=DIR)
rule quant:
        output:
                quant_file="{wd}/{sample}/salmon_out/quant.sf"
        run:
                #Path to quant file
                outfile=str(output.quant_file)
                #get srrid{sample}
                srrid=outfile.split("/")[1]
                #Run Salmon on sra object and delete fastq when finished.
                sra.SRA(srrid,directory=DIR).trim(trim_galore).quant(salmon).delete_fastq()
                
rule merge:
        input:
                ["{wd}/{sample}/salmon_out/quant.sf".format(wd=DIR,sample=s) for s in ra_ids]
        output:
                "{wd}/results_TPM_tx.tsv"
        run:
                #read in 1st quant file
                with open(input[0]) as f:
                        thisdata=f.read().splitlines()
                #Remove header
                thisdata.pop(0)
                #Run accession IDs
                names=[]
                #Transcript and gene IDs
                txids=[]
                gnids=[]

                #Get transcript and gene ids
                for l in thisdata:
                    thistx=l.split('\t')[0]
                    txids.append(thistx.split('|')[0])
                    gnids.append(thistx.split('|')[3])

                dftpm=pd.DataFrame({'Rat_ensembl_transcript_id_version':txids,'Rat_ensembl_gene_id_version':gnids})
                dfcount=pd.DataFrame({'Rat_ensembl_transcript_id_version':txids,'Rat_ensembl_gene_id_version':gnids})


                #Loop through all quant files
                for qf in input:
                        # current RAid
                        name=qf.split('/')[1]
                        names.append(name)
                        #Get TPM
                        tpm=pd.read_csv(qf,sep='\t',usecols=[3],skiprows=0)
                        #Get Counts
                        counts=pd.read_csv(qf,sep='\t',usecols=[4],skiprows=0)
                        #Add TPM and counts to respective df with current Run accession ID as column name
                        dftpm[name]=tpm['TPM']
                        dfcount[name]=counts['NumReads']

               
                
    
                ##transcript level TPMs and counts.
                dftpm.to_csv(wildcards.wd+"/results_TPM_tx.tsv",sep='\t',index=False)
                dfcount.to_csv(wildcards.wd+"/results_Counts_tx.tsv",sep='\t',index=False)


                #Collapse so that each gene id is listed once. sum up corresponding transcript TPM and counts.
                dftpm.drop('Rat_ensembl_transcript_id_version', axis=1, inplace=True)
                dfcount.drop('Rat_ensembl_transcript_id_version', axis=1, inplace=True)               
                dftpm = dftpm.groupby(['Rat_ensembl_gene_id_version'],as_index = False).sum()
                dfcount = dfcount.groupby(['Rat_ensembl_gene_id_version'],as_index = False).sum()



                dftpm.to_csv(wildcards.wd+'/results_TPM_gene.tsv',sep='\t',index=False)
                dfcount.to_csv(wildcards.wd+'/results_Count_gene.tsv',sep='\t',index=False)
