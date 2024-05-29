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
            "input/processed_sequences/16s/tables/full_tax_table_16s.txt",
            sep="\t", quote = FALSE, col.names=NA)

write.table(t(otu_table(ps)),
            "input/processed_sequences/16s/tables/full_seq_table_16s.txt",
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



