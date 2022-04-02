library(biomaRt)  
library(DESeq2)
library(dplyr)



#Read in count information.
countData = read.table("out/results_Count_gene.tsv",header=TRUE,sep = '\t',row.names="ensembl_gene_id_version")



#Round to nearest int
countData=round(countData,0) 


##Read in expermental design
metadata = read.table("Design.tsv",header=TRUE,row.names=1,sep = '\t')

#Should return TRUE
#all(rownames(metadata) == colnames(countData))


##Make DEseq2 object
dds = DESeqDataSetFromMatrix(countData = countData,colData = metadata,design = ~ Sample)
dds = DESeq(dds)
#Contrast case vs control
result = results(dds, contrast=c("Sample","case","control"))
## Remove rows with NA
result = result[complete.cases(result),]
#Put GeneID as column
result = cbind(ensembl_gene_id_version = rownames(result), result)




################################################################################

#Set up marts
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
rat = useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")


#Read in ensembl_transcript_id_version IDs
ID_List = read.table("out/results_TPM_tx.tsv",header=TRUE)
ID_List=as.vector(unlist(ID_List['ensembl_transcript_id_version']))

#Gather Rat data
df_Rat=getBM(attributes = c("ensembl_transcript_id_version","ensembl_gene_id_version","ensembl_gene_id","gene_biotype","description"), 
                 filters = "ensembl_transcript_id_version", 
                 values = ID_List , 
                 mart = rat)
                 
#Gather RatHuman Ortholog data                 
df_RatHumanOrtho = getLDS(attributes = c("ensembl_transcript_id_version","hsapiens_homolog_orthology_type"), 
                 filters = "ensembl_transcript_id_version", 
                 values = ID_List , 
                 mart = rat, 
                 attributesL = c("hgnc_symbol","ensembl_gene_id","gene_biotype","description"),             
                 martL = human, 
                 uniqueRows=T)

#Rename so that dfs have matching ensembl_transcript_id_version column                   
names(df_RatHumanOrtho)[names(df_RatHumanOrtho) == "Transcript.stable.ID.version"] <- "ensembl_transcript_id_version"               
#Merge dfs by ensembl_transcript_id_version
RatHumanOrthoMerged= merge(df_RatHumanOrtho,df_Rat,by="ensembl_transcript_id_version",all=TRUE)


##Drop transcript columns, unique rows, and merge.
RatHumanOrthoMerged=RatHumanOrthoMerged %>% select(-contains("ensembl_transcript_id_version"))
RatHumanOrthoMerged=distinct(RatHumanOrthoMerged)

#Merge with DESeq2 
RatHumanOrthoDGE = merge(RatHumanOrthoMerged,result,by="ensembl_gene_id_version",all.y=TRUE)
write.table(RatHumanOrthoDGE,"RatHumanOrthoDGE.tsv" ,sep = '\t',row.names = FALSE)

#Merge with counts. 
#Put GeneID as column
countData = cbind(ensembl_gene_id_version = rownames(countData), countData)
RatHumanOrthoCounts = merge(RatHumanOrthoMerged,countData,by="ensembl_gene_id_version",all=TRUE)
write.table(RatHumanOrthoCounts,"RatHumanOrthoCounts.tsv" ,sep = '\t',row.names = FALSE)


