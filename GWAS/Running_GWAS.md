# GWAS (ANGSD) and permutations

This pipeline outlines the necessary steps to run the GWAS and permutations.
The code below is not runnable in its current state, as the size of VCF files will require running the pipeline on a computing cluster, which vary widely in specifications.

## Input

VCFs for 479 individuals from Meier et al 2020

For each species there are two important files (which can be found in ./1.data/association/) , which have as many rows as individuals in the sample: 

* cov.txt - covariates for association (wing area, sex code (0=male, 1=female), admixture proportion) 
* cornea.area.txt - aspect ratio (elongatedness)

Make sure VCF order of individuals matches covariates order (I always use ascending order of unit/cam IDs). Check order of vcf with bcftools query -l x.vcf

Remove missing individuals (no phenotypes): We have 479 of 484.
```bash
bcftools view -S ^missing_covariate_ind.txt -Oz -o erato_filtered.recode.vcf.gz erato.vcf.gz

#Index
tabix -p vcf erato_filtered.recode.vcf.gz

#Check scaff length
bcftools query -f '%CHROM\n' erato_filtered.recode.vcf.gz | sort | uniq

```

## 1. NGSadmix

Already done by Montejo-Kovacevich et al. (2021). Obtain admixture proportions to control for population structure in the GWAS. Install NGSadmix from [here](http://www.popgen.dk/software/index.php/NgsAdmix), 

```bash
module load htslib/1.2.1 bcftools-1.9-gcc-5.4.0-b2hdt5n samtools/1.3.1 vcftools

# by Joana Meier: 
# replace with your file prefix
file=PC062_merged.PL.AD.HAPCUT2.inf0.05.corr.LDpruned
# Fix the missing genotypes (they need a PL)
zcat $file.vcf.gz | sed -e 's|./.:.:.:.:.:.|0/0:./.:0.33,0.33,0.34:0:10,10,10:0,0|g' | gzip > $file.corr.vcf.gz

# Generate a Beagle file
angsd -out $file -nThreads 8 \
  -vcf-gl $file.corr.vcf.gz -doGlf 2 \
  -doMajorMinor 1 -doMaf 2 -SNP_pval 1e-6 

# Run NGSadmix
nice NGSadmix -likes $file.beagle.gz -K 2 -minMaf 0.05 -outfiles $file.maf0.05.K2 

# Get 1 column for GWAS, which is added to cov.txt file
cut -f 1 -d " " $file.maf0.05.K2.qopt > NGSadmix.prop

```

## 2. ANGSD association

Current version of ANGSD has fixed vcf reading [angsd version: 0.933-101-g7ea6e4b (htslib: 1.10.2-131-g0456cec) build(Aug 24 2020 10:35:22)]- however, you have to use the option to **-vcf-pl** to specify we have PL likelihoods in our vcfs.

```bash
module load python-3.5.2-gcc-5.4.0-rdp6q5l gcc
module load R

REF=~/Heliconius_erato_demophoon_v1_-_scaffolds.fa
FILE=erato

angsd -yQuant ./1.data/assocation/era.479/cornea.area.txt \
 -doAsso 6 -doPost 1 -out out/${FILE}.gwas.shape.cov \
 -doMajorMinor 4 -doMaf 1 -vcf-pl ./erato.vcf.gz \
 -P 8 -model 1 -cov ./1.data/assocation/era.479/cov.txt \
 -ref $REF

# Compute window median, min p-values (see short explanation below)
Rscript ./3.scripts/extra/compute50SNP10SNP.windows.R out/${FILE}.gwas.shape.cov.lrt0.gz

# output will be "erato.gwas.shape.cov.lrt0.gz.50snp.10snp.pos"

```

The script 3.scripts/extra/compute50SNP10SNP.windows.R takes the output from ANGSD, calculates p.values from likelihood ratio test statistic (chi2 distribution, 1df) and creates sliding windows of 50SNP and 10SNP steps (fixed) with the R package winscanr ([](https://github.com/tavareshugo/WindowScanR)[tavareshugogithub](https://github.com/tavareshugo)), to obtain min, max, and median values from two columns: LRT.pval (p.values per snp) and position (to obtain mininum and maximum position)


## 3. Genome-wide permutations

Create 200 random “phenotype files” (i.e. cornea.area.txt) by randomly permuting the observed aspect ratio values. I use R:

```R
era.cov <-  read.csv("~/1.data/assocation/era.479/era.n479.covariates.csv")

era.ar <- era.cov$log.cornea.area

era.rand.ar <-list()
for(i in 1:200){
  era.rand.ar[[i]] <- sample(era.ar) }

for(i in 1:200){
  write.table(era.rand.ar[[i]], paste0("permutations/era.n479/cornea.area.", str_pad(i, 3, pad = "0"),".txt" ), row.names = F,col.names = F ) }

```

This will give rise to 200 files per species, called “permutations/era.n479/cornea.area.[001-200].txt”. These will be used as input for ANGSD.

On a cluster. Array job running GWAS 200 times per species, one for each permuted phenotype file
R script computed window averages for 50SNP win 10SNP step (same as observed data)

```bash
# create script perm.assoc.sh

REF=~/Herato_final_corrected.fa
PERM=`sed -n -e "$SLURM_ARRAY_TASK_ID p" permutations/era.n479/cornea.area.`

angsd -yQuant permutations/era.n479/cornea.area.${PERM}.txt \
-doAsso 6 -doPost 1 -out out/perm_${PERM}.gwas.shape.cov \
-doMajorMinor 4 -doMaf 1 -vcf-pl ~/erato_filtered.recode.vcf.gz \
-P 8 -model 1 -cov ~/1.data/assocation/era.479/cov.txt \
-ref $REF -fai ${REF}.fai -nInd 484

# Compute window medians, min p-values, same as with observed data
Rscript ~/3.scripts/extra/compute50SNP10SNP.windows.R out/perm_${PERM}.gwas.shape.cov.lrt0.gz

# to run array job
sbatch --array 1-200 perm.assoc.sh

```

### Inspect results

```
zcat erato_filtered.recode.gwas.shape.cov.lrt0.gz | head -n 20

#check all scaffs are present
zcat erato_filtered.recode.gwas.shape.cov.lrt0.gz | awk '{print $1}' | sort | uniq


head erato_filtered.recode.gwas.shape.cov.lrt0.gz.50snp.10snp.pos


```



### Rank observed p-values among empirical (permutated) p-values

First extract all median p-values of all windows and permutations.

Then the R script 3.scripts/extra/perm.ranking.R will import the median p-values from 200 permutations (200 values per outlier window), add the observed values (real p-values), and rank all p-values per window (observed +permuted, total=201). The observed value ranking will determine where in the null distribution it falls, we will only consider significant those above the 99th percentile (i.e. with a ranking of 1/201 or 2/201). 


```bash
# extract the column (i.e. all windows/rows) with median p-values from each permutation (careful col names may be shifted)`
for i in perm_*.gwas.shape.cov.lrt0.gz.50snp.10snp.pos
do
n=${i:17:3}
awk -F '\t' '{ print $9}' $i | column -t > lrt.median.all.rows.perm$n.txt
done

# concatenate all rows of all permutations
paste -d '\t'  lrt.median.all.rows.perm*.txt > cov3.all.rows.lrtperm.txt

# check correct number of columns, 200
awk -F'\t' '{print NF; exit}' cov3.all.rows.lrtperm.txt

# run r in a cluster with:
# input 1- concatenated permutation results and
# input 2- observed gwas results
nice Rscript ../../../../scripts/perm.ranking.R  cov3.all.rows.lrtperm.txt erato.gwas.shape.cov.lrt0.gz.50snp.10snp.pos &> perm.ranking.out&


```

### Plot the candidates in Manhattan plot

```R
library(dplyr)
library(stringr)
library(ggplot2)
 
rank_file <- "~/GWAS_project/erato_filtered.recode.gwas.shape.cov.lrt0.gz.50snp.10snp.pos.ranks"
 
gwas_ranks <- read.table(rank_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
 
hist(gwas_ranks$outlier.ranking.permutation200)

#create a unique SNP ID
gwas_ranks$SNP<-paste("r",1:length(gwas_ranks$scaff), sep="")

#Get position from scaff old position

# Read scaff - size file
chr_size_file="Herato_final_corrected.fa.chr_size.txt"
chr_size=read.table(chr_size_file, col.names=c("chr", "size"))
chr_size=chr_size %>% mutate(cum_sum = cumsum(size))

# format df
dxy.tmp <- gwas_ranks %>% 
  
  # Compute scaff size
  dplyr::group_by(scaff) %>% 
  dplyr::summarise(chr_len=max(Position_min)) %>% 
  
  # Calculate cumulative position of each scaff
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  dplyr::left_join(gwas_ranks, ., by=c("scaff"="scaff")) %>%
  
  # Add a cumulative position of each SNP
  arrange(scaff, Position_min) %>%
  mutate(BPcum=Position_min+tot) 

# get scaff center positions for x-axis
axisdf <- dxy.tmp %>% dplyr::group_by(scaff) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

# First, had to obtain it per scaffold
axisdf=split(dxy.tmp, dxy.tmp$scaff)
axisdf=lapply(axisdf, function(x) (max(x$BPcum) + min(x$BPcum))/2)
axisdf=as.data.frame(t(as.data.frame(bind_rows(axisdf))))
names(axisdf)="center"
axisdf$scaffold=row.names(axisdf)
axisdf$chr=str_sub(axisdf$scaffold, end=-3)

# Now per scaff (I hope)
axisdf=split(axisdf, axisdf$chr)
axisdf=lapply(axisdf, function(x) (max(x$center) + min(x$center))/2)
axisdf=as.data.frame(t(as.data.frame(bind_rows(axisdf))))
names(axisdf)="center"
axisdf$label=row.names(axisdf)

# Remove last two characters from each scaffolds, to gather those belonging to the same scaff
chr=unique(gwas_ranks$scaff)
chr=str_sub(chr, end=-3)

# Colors for plot, to identify scaffs
mypalette=c("blue4", "orange3")
mypalette=c("grey75", "grey30")

# Remove last two characters from each scaffolds, to gather those belonging to the same scaff
chr=unique(gwas_ranks$scaff)
chr=str_sub(chr, end=-3)

# List unique scaffs
chr_unique=unique(chr)
# Create colors vector by repeating the color sequence in mypalette
cols_seq=rep(mypalette, length(chr_unique)/length(mypalette))
# Bind chr name and color into a single df
chr_col=as.data.frame(cbind(chr_unique, cols_seq))
names(chr_col)=c("chr", "cols_seq")

# Turn chr into a data frame
chr=as.data.frame(chr)

# Merge df to associate chr with same color at each repeated chr value
chr_col=merge(chr, chr_col, by="chr")
ylims=c(0,0.5)

#Get outliers from full dataframe
outliers <- subset(dxy.tmp, dxy.tmp$outlier.ranking.permutation200 == 1 | dxy.tmp$outlier.ranking.permutation200 == 2)

#outliers <- subset(dxy.tmp, -log10(dxy.tmp$LRT.pval_median)>0.6)

dim(outliers)
[1] 849 17  #849 genes have a rank of 1 or 2. 

#make a list of the candidates
candidates <- dxy.tmp[dxy.tmp$SNP %in% outliers$SNP, ]$SNP

# Add a new column `Match` to df1 using base R
dxy.tmp$candidate <- ifelse(dxy.tmp$SNP %in% outliers$SNP, "candidate", "not_candidate")

#Plot
#refer to that list of candidates with the highlight parameter

p <- ggplot()+
  # Show all points
  geom_point(data=dxy.tmp, aes(x=BPcum, y=-log10(LRT.pval_median), color=as.factor(scaff)), alpha=0.5, size=0.5)+
  scale_color_manual(values = chr_col$cols_seq) +
  # custom X axis:
  scale_x_continuous( label = gsub("Herato", "", axisdf$label), breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0))+ #, limits = ylims) + # expand=c(0,0)removes space between plot area and x axis
  labs(x = NULL, y = "-log10(p)") +
  # Custom the theme:
  theme_classic() +
  theme(legend.position="none",
  axis.text = element_text(size = 12),
  axis.title = element_text(size = 12),
        axis.line.x = element_line(color="gray48", size = 0.4),
        axis.line.y = element_line(color="gray48", size = 0.4),
        plot.margin = unit(c(1, 1, 1, 1), "cm")) + ylim(c(0,1.2))

p <- p + 
  geom_point(data = outliers, aes(x = BPcum, y = -log10(LRT.pval_median)), color = "green", size=0.6) + 
  scale_x_continuous(labels = gsub("Herato", "", axisdf$label), breaks = axisdf$center) +
  scale_y_continuous(expand = c(0, 0)) + # Expand removes space between plot area and x-axis
  labs(x = NULL) +
  # Customize the theme
  theme_classic() +
  theme(legend.position = "none",
        axis.line.x = element_line(color = "gray48", size = 0.4),
        axis.line.y = element_line(color = "gray48", size = 0.4)) + ylim(c(0,1.2)) + geom_hline(yintercept=0.6, color = "blue")

#save plot
# Save the plot as a PDF
ggsave("manhattan_with_candidates_above06.pdf", plot = p, width = 7, height = 4)

#Or we also use the manhattan function...

manhattan(dxy.tmp, chr ="scaff", bp = "BPcum", p = "LRT.pval_median", snp = "SNP", highlight = candidates,
          col = c("blue4", "orange3"), main = "Manhattan Plot", xlab = "scaffold")
          
```


### For SI: extract lowest p-value per permutation

Grab minimum p-value per permutation (genome-wide), store. The 95th percentile will give a genome-wide threshold at p=0.05, presented in the SI Fig. S11.

Not very realistic for per site null distributions, as in every gwas of 25 million snps there will be some spuriously low p-values in each genome-wide permutation. 

```bash
# use p value min per window, column8, equivalent to finding minimum p-value genome wide (but faster)
for i in perm_*.gwas.shape.cov.lrt0.gz.50snp.10snp.pos; do awk 'NR==1 || $8 < min {line = $0; min = $8} END {print min}' $i >> min.LRT.pval.window.txt; done &

percentile_value=$(sort -n pvalues.txt | awk 'BEGIN {p=0.95} {a[NR]=$0} END {print a[int(NR*p)]}')
echo "The 95th percentile value is: $percentile_value"
#The 95th percentile value is: 4.99167985176457e-08

```

### Getting functional annotation from reference genome

```
 grep "evm.model.Herato1001.120" Heliconius_erato_demophoon_v1_gene_name_only_model_sub_sort_Herato1001.fa.tsv | cut -c 1-200

```










