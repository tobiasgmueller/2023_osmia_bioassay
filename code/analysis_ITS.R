# this is a scrip to follow the dada2_pipine_16s.R pipeline


#read in packages ####

if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
if(!requireNamespace("phyloseq")){
  install.packages("phyloseq")
}
if(!requireNamespace("DESeq2")){
  install.packages("DESeq2")
}



library(phyloseq); packageVersion("phyloseq")
library(tidyverse)
library(DESeq2)


# read in data from pipeline
rm(list=ls())


count_tab <- read.csv("input/its/asv_counts_ITS.csv", header=T, row.names=1)[ , -c(1:4)]

tax_tab <- as.matrix(read.csv("input/its/asv_its_taxonomy_nocontam.csv", header=T,
           row.names=1))

sample_info_tab <- read.csv("input/sequence_metadata.csv", row.names=1) %>% 
      mutate(across( c(treatment, day, description), as.factor) ) 
  



# normalize for sampling depth using a variance stabilizing transformation



# first we need to make a DESeq2 object
deseq_counts <- DESeqDataSetFromMatrix(count_tab, colData = sample_info_tab, design = ~type) 


deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)

    # and here is pulling out our transformed table
vst_trans_count_tab <- assay(deseq_counts_vst)

    # and calculating our Euclidean distance matrix
euc_dist <- dist(t(vst_trans_count_tab))











asv_tab <- read.csv("input/its/asv_counts_ITS.csv")

#set column 1 to row names
asv_tab<- asv_tab %>%
  `row.names<-`(., NULL) %>% 
   column_to_rownames(var = "X")


asv_tax <- read.csv("input/its/asv_taxonomy_ITS.csv")

#set column 1 to row names
asv_tax<- asv_tax %>%
  `row.names<-`(., NULL) %>% 
   column_to_rownames(var = "X")


asv_fasta <- readDNAStringSet("input/its/asv_ITS.fa")


#set column 1 to row names
asv_tax<- asv_tax %>%
  `row.names<-`(., NULL) %>% 
   column_to_rownames(var = "X")













# read in dada2 16s pipeline output 

taxa <- readRDS("input/processed_sequences/ITS/taxa_ITS.rds")
seqtab.nochim <- readRDS("input/processed_sequences/ITS/seqtab_nochim.rds")



# prep df about samples with id as the row names
samdf <- read_csv("sequencing_results/sequence_metadata.csv")%>%
  mutate(across(where(is.character), as.factor),
         across(where(is.numeric), as.factor)) %>% 
  column_to_rownames(var = "id")
  



# then construct phyloseq object

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

tax_table(ps) <- cbind(tax_table(ps),
                           rownames(tax_table(ps)))


# adding unique ASV names
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
tax_table(ps) <- cbind(tax_table(ps),
                           rownames(tax_table(ps)))
head(taxa_names(ps))


colnames(tax_table(ps)) <- c("Kingdom", "Phylum", "Class", "Order",
    "Family", "Genus", "Species", "ASV_SEQ", "ASV_ID")
ps


# write these raw tables

write.table(tax_table(ps),
            "input/processed_sequences/ITS/full_tax_table_ITS.txt",
            sep="\t", quote = FALSE, col.names=NA)

write.table(t(otu_table(ps)),
            "input/processed_sequences/ITS/full_seq_table_ITS.txt",
            sep="\t", quote = FALSE, col.names=NA)


# then save the phyloseq object
saveRDS(ps, "input/processed_sequences/ITS/ps_ITS.rds")








