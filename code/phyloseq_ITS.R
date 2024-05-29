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

taxa <- readRDS("input/processed_sequences/ITS/ITS_taxa.rds")
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





# plot richness
plot_richness(ps, x="treatment", measures=c("Shannon", "Simpson"), color="day")



# make an NMDS

ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")


plot_ordination(ps.prop, ord.nmds.bray, color="day", title="Bray NMDS")




# make bar chart

top100 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:100]
ps.top100 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top100 <- prune_taxa(top100, ps.top100)
plot_bar(ps.top20, x="treatment", fill="Class") + facet_wrap(~day, scales="free_x")


plot_bar(ps.top20, x="treatment", fill="Genus") 








