library("biomaRt")  

#Any quant.sf file
quantfile = read.table("out/SRR9164625/salmon_out/quant.sf",header=TRUE)
IDList=as.vector(unlist(PanHuman['Name']))





human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
rat = useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")


#Gather Rat data
df_Rat=getBM(attributes = c("ensembl_transcript_id_version","ensembl_gene_id","gene_biotype","description"), 
                 filters = "ensembl_transcript_id_version", 
                 values = GeneList , 
                 mart = rat)
                 
#Gather RatHuman Ortholog data                 
df_RatHumanOrtho = getLDS(attributes = c("ensembl_transcript_id_version","ensembl_gene_id","hsapiens_homolog_orthology_type"), 
                 filters = "ensembl_transcript_id_version", 
                 values = GeneList , 
                 mart = rat, 
                 attributesL = c("hgnc_symbol","ensembl_gene_id","gene_biotype","description"),             
                 martL = human, 
                 uniqueRows=T)
                 
names(df_RatHumanOrtho)[names(df_RatHumanOrtho) == "Transcript.stable.ID.version"] <- "ensembl_transcript_id_version"               

metadata= merge(df_RatHumanOrtho,df_Rat,by="ensembl_transcript_id_version",all=TRUE)

write.table(RatHumanOrtho,"RatHumanOrtho.tsv",sep = "\t")





                 
