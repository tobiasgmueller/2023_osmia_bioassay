# script for ITS dada 2 pipeline
#following DADA2 tutorial at https://benjjneb.github.io/dada2/ITS_workflow.html



#clear workspace
rm(list = ls())


#install code for dada 2
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")

if (!requireNamespace("dada2", quietly=TRUE))
  BiocManager::install("dada2")

if (!requireNamespace("DESeq2", quietly=TRUE))
  BiocManager::install("DESeq2")

if (!requireNamespace("decontam", quietly=TRUE))
  BiocManager::install("decontam")



#load packages and functions
library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
library(decontam)


# set path to fastq files
path <- "sequencing_results/ITS"  
list.files(path)

#split into forward and reverse
fnFs <- sort(list.files(path, pattern = "_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fastq.gz", full.names = TRUE))

# read in primer sequences to check for primers
FWD <- "CTTGGTCATTTAGAGGAAGTAA"  ## ITS1 forward
REV <- "GCTGCGTTCTTCATCGATGC"  ## ITS2R reverse

# function to find complements of DNA strings
allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
        RevComp = Biostrings::reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
}

# find all complements of primers
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

# filter out files with ambiguous bases before checking for primers
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))

filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)




# check how many times the primers appear 
primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), FWD.ReverseReads = sapply(FWD.orients,
    primerHits, fn = fnRs.filtN[[1]]), REV.ForwardReads = sapply(REV.orients, primerHits,
    fn = fnFs.filtN[[1]]), REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))




# cut out remaining primers ####
# there are a few primers. lets remove these using Cutadapt - which for windows requires a C++ compiler

#if running in windows set to cutadapt path
#cutadapt <- "C:/path/cutadapt" # CHANGE ME to the cutadapt path on your machine
#system2(cutadapt, args = "--version") # Run shell commands from R
#then update the below system2 calls to get rid of quotes on cutadapt

system2("cutadapt", args = "--version") # Run shell commands from R


# and now trim primers
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2("cutadapt", args = c(R1.flags, R2.flags,
                             "-m", 20, # drop reads shorter than 20 to remove 0 length reads
                             "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

#check that cutting worked
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), FWD.ReverseReads = sapply(FWD.orients,
    primerHits, fn = fnRs.cut[[1]]), REV.ForwardReads = sapply(REV.orients, primerHits,
    fn = fnFs.cut[[1]]), REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))



# now that primers are cut, we pair forward and reverse reads
cutFs <- sort(list.files(path.cut, pattern = "_R1.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2.fastq.gz", full.names = TRUE))


# Extract sample names
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)


# inspect read quality
plotQualityProfile(cutFs[1:3])
plotQualityProfile(cutRs[1:3])

#plotQualityProfile(cutFs, aggregate = TRUE)
#plotQualityProfile(cutRs, aggregate = TRUE)

ShortRead::readFastq(cutFs[1])
ShortRead::readFastq(cutRs[1])


# filter and trim ####

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))


out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs,
              maxN=0,
              maxEE=c(2,2),
              truncQ=2,
              minLen = 50,
              rm.phix=TRUE,
              compress=TRUE,
              multithread=TRUE) # On Windows set multithread=FALSE
head(out)


# learn error rate ####
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#then plot error rate
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# sample inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE, pool="pseudo")
dadaRs <- dada(filtRs, err=errR, multithread=TRUE, pool="pseudo")


# inspect that it went okay
dadaFs[[1]]
dadaRs[[1]]


#Merge paird reads ####
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])




# construct ASV table ####
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

table(nchar(getSequences(seqtab)))




# remove chimeras ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

# percent of merged sequence reads that are chimeras
sum(seqtab.nochim)/sum(seqtab)




# track reads through pipeline ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#save track
write.csv(track, file="sequencing_results/ITS/track_through_pipe.csv")


# assign taxonomy using UNITE database
unite.ref <- "sequencing_results/ITS/tax/sh_general_release_dynamic_04.04.2024.fasta"  
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)


# look at taxanomic assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)






#then for save "raw" outputs
saveRDS(taxa, "input/ITS/taxa_ITS.rds")
saveRDS(seqtab.nochim, "input/ITS/seqtab_nochim_ITS.rds")


asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
    asv_headers[i] <- paste("ASV", i, sep="_")
}


# fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "input/ITS/asv_ITS.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- asv_headers
write.csv(asv_tab, "input/ITS/asv_counts_ITS.csv")

#taxa table
asv_taxa<-taxa
row.names(asv_taxa) <- asv_headers
write.csv(asv_taxa, "input/ITS/asv_taxonomy_ITS.csv")



# now check for contaminants

# create vector saying which samples are controls
vector_for_decontam <- c(rep(TRUE, 1), rep(FALSE, 30))

contam_df <- isContaminant(t(asv_tab), neg=vector_for_decontam)

table(contam_df$contaminant) # identified 6 as contaminants

# getting vector holding the identified contaminant IDs
contam_asvs <- row.names(contam_df[contam_df$contaminant == TRUE, ])

# in this case its an unidentified fungi
asv_tax[row.names(asv_tax) %in% contam_asvs, ]



  

# write out decontaminated files


  # making new fasta file
contam_indices <- which(asv_fasta %in% paste0(">", contam_asvs))
dont_want <- sort(c(contam_indices, contam_indices + 1))
asv_fasta_no_contam <- asv_fasta[- dont_want]

    # making new count table
asv_tab_no_contam <- asv_tab[!row.names(asv_tab) %in% contam_asvs, ]

    # making new taxonomy table
asv_tax_no_contam <- asv_tax[!row.names(asv_tax) %in% contam_asvs, ]

    ## and now writing them out to files
write(asv_fasta_no_contam, "input/its/asv_its_nocontam.fa")
write.csv(asv_tab_no_contam, "input/its/asv_its_counts_nocontam.csv")
write.csv(asv_tax_no_contam, "input/its/asv_its_taxonomy_nocontam.csv")







