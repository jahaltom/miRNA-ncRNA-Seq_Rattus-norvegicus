# BioProject PRJNA545590 (ncRNA-Seq)

# miRNA-Seq 
Since the data is paired-end, it is very unlikely to find mature miRNAs, the best thing that can happen is finding long precursors. Since miRDeep2 only does single end, I used the 1st mate of the pair which is standard.

Taking run accession IDs as input (RAids.txt) this pipeline quantifies miRNA-Seq data using miRDeep2.

The quantification step in this pipeline is set up to allow for parallelization with each run accession ID being done separately. Since a mature miRNA can come from more than 1 precursor miRNA, these particular mature miRNA's are averaged across precursor miRNA's. 

Species is Rattus norvegicus.


**Create and activate conda enviroment:**
```
conda env create -f environment.yaml
conda activate miRNA-Seq
```

**Download required fasta files from miRBase and unzip**
```
wget https://www.mirbase.org/ftp/CURRENT/mature.fa.gz
wget https://www.mirbase.org/ftp/CURRENT/hairpin.fa.gz
gunzip *gz*
```


**Run miRNA-Seq quantification analysis**
This should generate "countTable.tsv" which is used as count input for DESeq2.

```
snakemake -j 6 -s miRNA-Seq.py --cluster "sbatch -t 02:00:00  -c 30 -p RM-shared"
```
Findings: Failed due to lack of data. DESeq2 not possible.
```
SRR9164621
#miRNA  read_count      precursor       total   621     621(norm)
rno-miR-466b-3p 1.00    rno-mir-466b-1  1.00    1.00    500000.00
rno-miR-466b-3p 1.00    rno-mir-466b-3  1.00    1.00    500000.00
```



# mRNA/ncRNA-Seq

**Prepareing Transcriptome Data**

Creates a salmon index for all transcripts of the Rat Ensembl genes, including ncRNA. Uses Rat genome as a decoy file. 
```
sbatch Prepare_Transcriptome_data.sh
```

**Transcript Quantification**
Quantifies runs in parallel. Outputs transcript and gene level TPM and counts. Gene level summed up across transcripts.
```
snakemake -j 50 -s SRA_Quant.py --cluster "sbatch -t 03:00:00 -c 16 -N 1"
```
**DESeq2**
```
Rscript DESeq2_RatHumanOrthologMetadata.r
```

