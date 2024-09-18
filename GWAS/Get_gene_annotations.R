##!Rscript --vanilla
# Load packages required

# Elisa 2024 Modified from https://onlinelibrary.wiley.com/doi/10.1111/mec.16067 
#packages

library(data.table)
library(GenomicRanges)
library(GenomicFeatures)
library(ggplot2)
library(dplyr)
library(tidyr)
library(qqman)
options(scipen = 999)
library(stringr)

#0. data prep
# use 10 snp no overlaps to look at densities
era.shape.gwas <- read.table("outliers_rank_1_2.txt", header = T)

# references available at LepBase
ref.scaff.era <- read.table("./Herato_final_corrected.fa.fai", row.names = NULL)

# 1. get outlier ranges 

# make granges, keep metadata, xxx blocks
head(era.shape.gwas)
era.all.summ.ranges <-  makeGRangesFromDataFrame(era.shape.gwas, keep.extra.columns = T, 
                                                 start.field="Position_min",end.field=c("Position_max"),
                         seqnames.field=c("scaff"),
                         ignore.strand =T)


#2. extract genes (all and outlier)
#get genes from gff with genomic ranges
# create txdb from gff https://www.biostars.org/p/167818/
txdb.era <- makeTxDbFromGFF("1.data/Heliconius_erato_demophoon_v1.gff3" , format="gff3")

# checks
exonsBy(txdb.era, by="gene")
mcols(GenomicFeatures::cds(txdb.era, columns=c("TXID", "TXNAME", "GENEID"))); GenomicFeatures::transcripts(txdb.era)
nrow(mcols(GenomicFeatures::genes(txdb.era, columns=c("TXID", "TXNAME", "GENEID")))) # total 13676 genes

# obtain gene lengths for all genes of each species, called "width"
era.gene.lengths <-as.data.frame(GenomicFeatures::transcripts(txdb.era))

#2.2 erato outlier genes
length(GenomicFeatures::genes(txdb.era)) 

# use txname to overlap with ranges, as this has model instead of TU
overlap.genes.era <- subsetByOverlaps( GenomicFeatures::genes(txdb.era, columns=c("TXNAME")), era.all.summ.ranges, type="any")

# Step 3: To preserve the LRT.pval data, you can find the overlaps manually and merge metadata:
# Find overlaps manually
hits <- findOverlaps(overlap.genes.era, era.all.summ.ranges)

# Extract metadata from era.all.summ.ranges based on hits
metadata <- mcols(era.all.summ.ranges)[subjectHits(hits),]

# Combine this metadata with the overlapping genes data
overlap.genes.era <- overlap.genes.era[queryHits(hits)]
mcols(overlap.genes.era) <- cbind(mcols(overlap.genes.era), metadata)

# Now overlap.genes.era has both the gene information and the LRT.pval data
head(overlap.genes.era)

# Step 4: Order by LRT.pval
ordered_genes <- overlap.genes.era[order(mcols(overlap.genes.era)$LRT.pval_median)]

# Now ordered_genes is sorted by LRT.pval
head(ordered_genes)

# Ensure unique names in GRanges object by adding a unique suffix if needed
unique_names <- make.unique(names(ordered_genes))
names(ordered_genes) <- unique_names

# Now convert to data.frame
ordered_genes_df <- as.data.frame(ordered_genes)

# Order by LRT.pval
ordered_genes_df <- ordered_genes_df[order(ordered_genes_df$LRT.pval_median), ]

# Function to convert list to a comma-separated string
convert_list_to_string <- function(x) {
  if (is.list(x)) {
    sapply(x, function(y) paste(unlist(y), collapse = ","))
  } else {
    x
  }
}

# Convert all list columns to character strings
ordered_genes_df <- data.frame(lapply(ordered_genes_df, convert_list_to_string), stringsAsFactors = FALSE)

# Print structure to verify conversion
str(ordered_genes_df)

# Save the cleaned data frame to a text file
write.table(ordered_genes_df, file = "ordered_genes_cleaned.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#####################################################################
head(as.data.frame(overlap.genes.era)); str(as.data.frame(overlap.genes.era$TXNAME))

# use values (as these will match real names) # 4k genes
genes.era <-  as.data.frame(overlap.genes.era$TXNAME)$value; length(genes.era); length(unique(genes.era)) #1081 genes, 624 unique genes !!!
# write out table
write.table(genes.era,file="subsetted.candidate.genes.txt", sep = '/t',quote=FALSE, col.name=FALSE)


# **NEW**: Read the candidate genes and merge with GWAS data
# Read the candidate genes into a dataframe
candidate_genes <- read.table("ordered.candidate.genes.txt", header=TRUE, col.names="GeneID")
genes_with_ID <- read.table("selected_genes.gff", header = F)

# Extract the relevant part from candidate genes (removing '1/' prefix and keeping the rest)
candidate_genes$ProcessedGeneID <- str_extract(candidate_genes$GeneID, "(?<=1/).+")
# Remove the first two characters from each ID
candidate_genes$ProcessedGeneID <- substr(candidate_genes$GeneID, 3, nchar(candidate_genes$GeneID))

# Extract the part of the ID after the first `;`
genes_with_ID$ProcessedGWASID <- str_extract(genes_with_ID$V9, "(?<=;)[^;]+")
# Remove the first 7 characters ("Parent=")
genes_with_ID$ProcessedGWASID <- str_sub(genes_with_ID$ProcessedGWASID, 8)

# Match based on the processed columns
merged_data <- merge(candidate_genes, genes_with_ID, by.x="ProcessedGeneID", by.y="ProcessedGWASID")

# Order by LRT.pval_median
ordered_results <- merged_data %>% arrange(LRT.pval_median)

# Write the ordered results to a file
write.table(ordered_results, file="ordered_candidate_genes_by_LRT.txt", sep='\t', quote=FALSE, row.names=FALSE)
