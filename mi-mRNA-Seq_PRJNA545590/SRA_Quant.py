import yaml
import sys
import os
from pyrpipe import sra,quant
import pandas as pd


#Read in config file
configfile: "config.yaml"
DIR = config['DIR']
THREADS=config['THREADS']
Tr=config['transcriptome']



#Read in run_accession ids
with open ("RAids.txt") as f:
        ra_ids=f.read().splitlines()

#create salmon object. If index is not present it will be created
salmon=quant.Salmon(index="human_data/salmon_index",transcriptome=Tr,threads=THREADS)

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
                sra.SRA(srrid,directory=DIR).quant(salmon).delete_fastq()

rule merge:
        input:
                ["{wd}/{sample}/salmon_out/quant.sf".format(wd=DIR,sample=s) for s in ra_ids]
        output:
                outfile="{wd}/results_TPM_tx.tsv"
        run:
                #read in quant file
                with open(input[0]) as f:
                        thisdata=f.read().splitlines()
                thisdata.pop(0)
                #Run accession IDs
                names=[]
                #Transcript and gene IDs
                txids=[]

                #Get transcript ids
                for l in thisdata:
                        thistx=l.split('\t')[0]
                        if '|' in thistx: thistx=thistx.split('|')[0]
                        txids.append(thistx)

                dftpm=pd.DataFrame({'TranscriptID':txids})
                dfcount=pd.DataFrame({'TranscriptID':txids})


                #current quant file
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

                ##transcript level TPMs.
                #Read in metadata
                tx_md=pd.read_csv('Transcript_level_metadata.tsv',sep='\t',skiprows=0)
                #Fetch transcript ids and TPMs for all RA ids.
                df_tx=dftpm[['TranscriptID']+names].copy()
                #Merge with metadata
                df_tx = tx_md.merge(df_tx, on=['TranscriptID'])
                df_tx.to_csv(output.outfile,sep='\t',index=False)


                #Read in gene metadata
                md=pd.read_csv('Gene_level_metadata.tsv',sep='\t',skiprows=0)
                df_gene_tpm=dftpm[['TranscriptID']+names].copy()
                df_gene_count=dfcount[['TranscriptID']+names].copy()

                #Make df with transcript id and corresponding Gene_stable_ID
                md2=tx_md[['TranscriptID','Gene_stable_ID']]
                ##Merge TPM and count data with ids
                df_gene_tpm=md2.merge(df_gene_tpm, on=['TranscriptID'], how='right')
                df_gene_count=md2.merge(df_gene_count, on=['TranscriptID'], how='right')


                #Collapse so that each gene id is listed once. sum up corresponding transcript TPM and counts.
                df_gene_tpm.drop('TranscriptID', axis=1, inplace=True)
                df_gene_count.drop('TranscriptID', axis=1, inplace=True)
                df_gene_tpm = df_gene_tpm.groupby(['Gene_stable_ID'],as_index = False).sum()
                df_gene_count = df_gene_count.groupby(['Gene_stable_ID'],as_index = False).sum()


                ##Merge metadata to counts and TPM
                df_gene_tpm=md.merge(df_gene_tpm, on=['Gene_stable_ID'], how='right')
                df_gene_count=md.merge(df_gene_count, on=['Gene_stable_ID'], how='right')

                #reorder so gene name is first.
                df_gene_tpm = df_gene_tpm[ ['Gene_name'] + [ col for col in df_gene_tpm.columns if col != 'Gene_name' ] ]
                df_gene_count = df_gene_count[ ['Gene_name'] + [ col for col in df_gene_count.columns if col != 'Gene_name' ] ]



                df_gene_tpm.to_csv(DIR+'/results_TPM_gene.tsv',sep='\t',index=False)
                df_gene_count.to_csv(DIR+'/results_Count_gene.tsv',sep='\t',index=False)
