library(biomaRt)  
library(DESeq2)
library(dplyr)



#Read in count information.
countData = read.table("out/results_Count_gene.tsv",header=TRUE,sep = '\t',row.names="ensembl_gene_id_version_Rat")



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
result = cbind(ensembl_gene_id_version_Rat = rownames(result), result)




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
#Add Rat Sufix to all columns
colnames(df_Rat) <- paste(colnames(df_Rat),"Rat",sep="_")   

              
#Gather RatHuman Ortholog data                 
df_RatHumanOrtho = getLDS(attributes = c("ensembl_transcript_id_version","hsapiens_homolog_orthology_type"), 
                 filters = "ensembl_transcript_id_version", 
                 values = ID_List , 
                 mart = rat, 
                 attributesL = c("hgnc_symbol","ensembl_gene_id","gene_biotype","description"),             
                 martL = human, 
                 uniqueRows=T)

#Add Human Sufix to all columns
colnames(df_RatHumanOrtho) <- paste(colnames(df_RatHumanOrtho),"Human",sep="_")
#Rename so that dfs have matching ensembl_transcript_id_version column.                   
names(df_RatHumanOrtho)[names(df_RatHumanOrtho) == "Transcript.stable.ID.version_Human"] <- "ensembl_transcript_id_version_Rat"    
names(df_RatHumanOrtho)[names(df_RatHumanOrtho) == "Human.homology.type_Human"] <- "Human.homology.type"            

#Merge dfs by ensembl_transcript_id_version
RatHumanOrthoMerged= merge(df_RatHumanOrtho,df_Rat,by="ensembl_transcript_id_version_Rat",all=TRUE)



##Drop transcript columns, unique rows, and merge.
RatHumanOrthoMerged=RatHumanOrthoMerged %>% select(-contains("ensembl_transcript_id_version_Rat"))
RatHumanOrthoMerged=distinct(RatHumanOrthoMerged)

#Merge with DESeq2 
RatHumanOrthoDGE = merge(RatHumanOrthoMerged,result,by="ensembl_gene_id_version_Rat",all.y=TRUE)
write.table(RatHumanOrthoDGE,"RatHumanOrthoDGE.tsv" ,sep = '\t',row.names = FALSE)

#Merge with counts. 
#Put GeneID as column
countData = cbind(ensembl_gene_id_version = rownames(countData), countData)
RatHumanOrthoCounts = merge(RatHumanOrthoMerged,countData,by="ensembl_gene_id_version_Rat",all=TRUE)
write.table(RatHumanOrthoCounts,"RatHumanOrthoCounts.tsv" ,sep = '\t',row.names = FALSE)
