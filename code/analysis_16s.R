# this is a scrip to follow the dada2_pipine_16s.R pipeline


#read in packages ####

if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}

if(!requireNamespace("phyloseq")){
  install.packages("phyloseq")
}

library(phyloseq); packageVersion("phyloseq")
library(tidyverse)
library(DESeq2)


# read in data from pipeline
rm(list=ls())


count_tab <- read.csv("input/16s/asv_counts_16s.csv", header=T, row.names=1)[ , -c(1:4)]

tax_tab <- as.matrix(read.csv("input/16s/asv_16s_taxonomy_nocontam.csv", header=T,
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

























# read in dada2 16s pipeline output 

taxa <- readRDS("input/processed_sequences/16s/taxa.rds")
seqtab.nochim <- readRDS("input/processed_sequences/16s/seqtab_nochim.rds")



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
            "input/processed_sequences/16s/full_tax_table_16s.txt",
            sep="\t", quote = FALSE, col.names=NA)


write.table(t(otu_table(ps)),
            "input/processed_sequences/16s/full_seq_table_16s.txt",
            sep="\t", quote = FALSE, col.names=NA)




# check for mitochondria, chloroplast, Eukaryotas

"Chloroplast" %in% tax_table(ps)
"Mitochondria" %in% tax_table(ps)
"Eukaryota" %in% tax_table(ps)

#this still has alot of mitochondrial and chloroplast reads so we'll remove those
# remove mitochondria and chloroplast



# generate a file with mitochondria ASVs
MT1 <- subset_taxa(ps, Family == "Mitochondria")
MT1 <-  as(tax_table(MT1), "matrix")
MT1 <- MT1[, (ncol(MT1))]
MT1df <- as.factor(MT1)
goodTaxa <- setdiff(taxa_names(ps), MT1df)
ps_no_mito <- prune_taxa(goodTaxa, ps)
ps_no_mito


CH1 <- subset_taxa(ps_no_mito, Order == "Chloroplast")
CH1 <-  as(tax_table(CH1), "matrix")
CH1 <- CH1[, (ncol(CH1))]
CH1df <- as.factor(CH1)
goodTaxa <- setdiff(taxa_names(ps_no_mito), CH1df)
ps_no_chloro <- prune_taxa(goodTaxa, ps_no_mito)
ps_no_chloro



# then set that to ps for easy coding
ps<-ps_no_chloro


# then save the phyloseq object
saveRDS(ps, "input/processed_sequences/16s/ps_16s.rds")

