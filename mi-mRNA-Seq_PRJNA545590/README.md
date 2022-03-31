# BioProject PRJNA545590 (ncRNA-Seq)

# miRNA-Seq 

Taking run accession IDs as input (RAids.txt) and an experimental design (Design.tsv), this pipeline quantifies miRNA-Seq data using miRDeep2 and then conducts differential expression analysis with DESeq2.

The quantification step in this pipeline is set up to allow for parallelization with each run accession ID being done separately. Since a mature miRNA can come from more than 1 precursor miRNA, these particular mature miRNA's are averaged across precursor miRNA's. The resulting average is then rounded to the nearest integer.

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
This generates "countTable.tsv" which is used as count input for DESeq2.

```
snakemake -j 6 -s miRNA-Seq.py --cluster "sbatch -t 02:00:00  -c 30 -p RM-shared"
```


**Run DE analysis**

Failed: Since the data is paired-end, it is very unlikely to find mature miRNAs, the best thing that can happen is finding long precursors. Since miRDeep2 only does single end, I used the 1st mate of the pair which is standard. Findings:



```
SRR9164621
#miRNA  read_count      precursor       total   621     621(norm)
rno-miR-466b-3p 1.00    rno-mir-466b-1  1.00    1.00    500000.00
rno-miR-466b-3p 1.00    rno-mir-466b-3  1.00    1.00    500000.00
```
