###PREPPING YOUR R ENVIRONMENT###

##Set your working directory##
#Use getwd() to find out where you currently are
getwd()
#Use setwd() to set your working directory 
setwd("C:/Users/bradshawd/Documents/Bioinformatic_Analysis/Clam_Feed/R")
#Rerun getwd() to check it worked
getwd()

#Copy for future use
clipr::write_clip(my_df)

##Load up SummarySE function##

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summarized
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

##Load packages##
library(phyloseq)
library(vegan)
library(ggplot2)
library(plyr)
library(tidyverse)
library(FSA)
library(dplyr)
library(reshape)
library(DESeq2)
library(metagenomeSeq)
library(ggvenn)

###GETTING YOUR DATA INTO R###

##Import files##

BIOM <- import_biom(file.choose()) #phyloseq.biom
TREE =  read_tree(file.choose()) #tree/tree.nwk
META <- import_qiime_sample_data(file.choose()) #clam_feed.txt
META2  = read.delim(file.choose(), row.names=1) #clam_feed.txt

#check that sample names are the same between metadata and abundance 
sample_names(META)
sample_names(BIOM)

#Merge three items into one phyloseq object called data#
data <- merge_phyloseq (BIOM,TREE,META)

###DATA MANIPULATION AND SUMMARIZING SEQUENCES### 

#See what present names are#
colnames(tax_table(data))

#Switch them to classic format#
colnames(tax_table(data))=c("Kingdom","Phylum","Class","Order","Family","Genus", "Species")

#Check that it worked#
colnames(tax_table(data))

#Check number of samples
nsamples(data) #48#

##Filtering ASVs##

#Check number of taxa#
ntaxa(data) #1576#

#Filter low abundance ASVs

Fdata = filter_taxa(data, function(x) sum(x) >10, TRUE)
ntaxa(Fdata) #988#

#Copy out asv abundance (phyloseq uses otu, probably because it is just an aesthetics thing) table as a separate matrix from phyloseq object into R environment#
otutable = as(otu_table(Fdata), "matrix")

#Bind to it the taxonomy table which you copy out from phyloseq in the same step via the same rownames#
taxotutable = cbind(as(otutable, "matrix"), as(tax_table(data)[rownames(otutable), ], "matrix"))

#Export ASV & Taxonomy table out of R using the write function. Can write it as many things, but csv is my preferred method, cause easier to open in Excel#
write.csv(taxotutable , "test.otu.tax.csv")

##Summarize sequences by ASVs 

#Unfiltered version first#

#export table of sequence sums#
data_seqs_per_ESV <- as.data.frame(taxa_sums(data))

#Change colnames to Sequences#
colnames(data_seqs_per_ESV) <- c("Sequences")
sum(data_seqs_per_ESV$Sequences) #187315#

#Use cbind to get taxonomy added in#
data_seqs_per_ESV = cbind(as(data_seqs_per_ESV, "data.frame"), as(tax_table(data)[rownames(data_seqs_per_ESV), ], "matrix"))

#export it as a csv, and tell it that the the rownames should be named OTUID#
write.csv(data.frame("OTUID" =rownames(data_seqs_per_ESV), data_seqs_per_ESV) , "seqs_per_ESV.csv", row.names=FALSE)

#Filtered version second#

Fdata_seqs_per_ESV <- as.data.frame(taxa_sums(Fdata))
colnames(Fdata_seqs_per_ESV) <- c("Sequences")
sum(Fdata_seqs_per_ESV$Sequences) #184217#
Fdata_seqs_per_ESV = cbind(as(Fdata_seqs_per_ESV, "data.frame"), as(tax_table(Fdata)[rownames(Fdata_seqs_per_ESV), ], "matrix"))
write.csv(data.frame("OTUID" =rownames(Fdata_seqs_per_ESV), Fdata_seqs_per_ESV) , "filtered_seqs_per_ESV.csv", row.names=FALSE)

##Summarize sequences by samples (filtered and unfiltered)##

##INFOMATION USED FOR SUPPLEMENTARY TABLE 1


#Unfiltered first
data_seqs_per_sample <- as.data.frame(sample_sums(data))
colnames(data_seqs_per_sample) <- c("Full_Sequences")
sum(data_seqs_per_sample$Full_Sequences) #187315#
Fdata_seqs_per_sample <- as.data.frame(sample_sums(Fdata))
colnames(Fdata_seqs_per_sample) <- c("Trimmed_Sequences")
sum(Fdata_seqs_per_sample$Trimmed_Sequences) #184217#

#Sum of sequences should be same as above. Instead of exporting two different tables you can instead combine them with cbind and export that#
Comdata_seqs_per_sample = cbind(as(data_seqs_per_sample, "data.frame"), as(Fdata_seqs_per_sample, "data.frame"))
write.csv(data.frame("OTUID" =rownames(Comdata_seqs_per_sample), Comdata_seqs_per_sample) , "sequences_per_sample.csv", row.names=FALSE)

##Determine numbers of various taxonomic levels: total and defined

##INFOMATION USED FOR SUPPLEMENTARY TABLE 2

#Total numbers are all the unique ids given to organisms at that level while defined means that the id is actually biologically relevant instead of being uncultured, or unknown, or a repeat of the previous taxonomic level

#Extract out the taxonomy table from filtered dataset
Fdatataxa <- as.data.frame(tax_table(Fdata))

#Remove the labels at each taxonomic level to allow comparisons between taxonomic levels without it being a compounding variable
#For examples f__Unknown Family and g__Unknown Family would not be comparable otherwise
Fdatataxa <- data.frame(lapply(Fdatataxa, function(x) {gsub("d__", "", x)}))
Fdatataxa <- data.frame(lapply(Fdatataxa, function(x) {gsub("p__", "", x)}))
Fdatataxa <- data.frame(lapply(Fdatataxa, function(x) {gsub("c__", "", x)}))
Fdatataxa <- data.frame(lapply(Fdatataxa, function(x) {gsub("o__", "", x)}))
Fdatataxa <- data.frame(lapply(Fdatataxa, function(x) {gsub("f__", "", x)}))
Fdatataxa <- data.frame(lapply(Fdatataxa, function(x) {gsub("g__", "", x)}))
Fdatataxa <- data.frame(lapply(Fdatataxa, function(x) {gsub("s__", "", x)}))

#Kingdom
#Total
Fdatataxa%>%distinct(Kingdom, .keep_all=TRUE)%>%nrow #2
#Defined
Fdatataxa%>%distinct(Kingdom, .keep_all=TRUE)%>% filter(!is.na(Kingdom), Kingdom!="uncultured") %>%nrow #2

#Phylum
#Total
Fdatataxa%>%distinct(Kingdom, Phylum, .keep_all=TRUE)%>%nrow #17
#Defined
Fdatataxa%>%distinct(Kingdom, Phylum, .keep_all=TRUE)%>% filter(!is.na(Phylum), Phylum!="uncultured") %>%nrow #17

#Class
#Total
Fdatataxa%>%distinct(Kingdom, Phylum, Class, .keep_all=TRUE)%>%nrow #25
#Defined
Fdatataxa%>%distinct(Kingdom, Phylum, Class, .keep_all=TRUE)%>% filter(!is.na(Class), !grepl(pattern = "uncultured", x = Class)) %>% filter(!grepl(pattern = "unidentified", x = Class))%>% filter(!grepl(pattern = "metagenome", x = Class))%>%filter(as.character(Phylum) != as.character(Class))%>%nrow #24

#Order
#Total
Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, .keep_all=TRUE)%>%nrow #66
#Defined
Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, .keep_all=TRUE)%>% filter(!is.na(Order), !grepl(pattern = "uncultured", x = Order)) %>% filter(!grepl(pattern = "unidentified", x = Order))%>% filter(!grepl(pattern = "metagenome", x = Order))%>%filter(!grepl(pattern = "Unknown", x = Order))%>%filter(as.character(Class) != as.character(Order))%>%nrow #65

#Family
#Total
Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, Family, .keep_all=TRUE)%>%nrow #116
#Defined
Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, Family, .keep_all=TRUE)%>% filter(!is.na(Family), !grepl(pattern = "uncultured", x = Family)) %>% filter(!grepl(pattern = "unidentified", x = Family))%>% filter(!grepl(pattern = "metagenome", x = Family))%>%filter(!grepl(pattern = "Unknown", x = Family))%>%filter(as.character(Order) != as.character(Family))%>%nrow #112

#Genus
#Total
Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, .keep_all=TRUE)%>%nrow #189
#Defined
Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, .keep_all=TRUE)%>% filter(!is.na(Genus), !grepl(pattern = "uncultured", x = Genus)) %>% filter(!grepl(pattern = "unidentified", x = Genus))%>% filter(!grepl(pattern = "metagenome", x = Genus))%>%filter(!grepl(pattern = "Unknown", x = Genus))%>%filter(as.character(Family) != as.character(Genus))%>%nrow #160

#Species
#Total
Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, Species, .keep_all=TRUE)%>%nrow #236
#Defined
Fdatataxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, Species, .keep_all=TRUE)%>% filter(!is.na(Species), !grepl(pattern = "uncultured", x = Species)) %>% filter(!grepl(pattern = "unidentified", x = Species))%>% filter(!grepl(pattern = "uncultured", x = Genus))%>% filter(!grepl(pattern = "metagenome", x = Species))%>%filter(!grepl(pattern = "Unknown", x = Species))%>%filter(as.character(Genus) != as.character(Species))%>% filter(!grepl(pattern = "bacterium", x = Species))%>% filter(!grepl(pattern = "enrichment", x = Species))%>% filter(!grepl(pattern = "sp.", x = Species))%>% filter(!grepl(pattern = "endosymbiont", x = Species))%>%nrow #47-7 blanks = 19

###ALPHA DIVERSITY###

##Calculate alpha diversity##

Falphadiv = estimate_richness(Fdata, split = TRUE)
write.csv(Falphadiv, file='Falphadiv.csv')
Falphadiv <- read.csv("Falphadiv.csv", row.names = 1)

#cbind in the metadata from the filtered phyloseq#
Falphadiv_metadata = cbind(as(Falphadiv, "data.frame"), as(sample_data(Fdata)[rownames(Falphadiv), ], "data.frame"))

##Test for normality##
shapiro.test(Falphadiv_metadata$Shannon)#0.2707
shapiro.test(Falphadiv_metadata$Observed)#0.03545
shapiro.test(Falphadiv_metadata$Fisher)#0.002568
shapiro.test(Falphadiv_metadata$Simpson)#0.7609

##Correlation Tests##
cor.test(Falphadiv_metadata$Shannon, Falphadiv_metadata$Observed, method = "spearman", exact = FALSE)#7.016e-16
cor.test(Falphadiv_metadata$Shannon, Falphadiv_metadata$Fisher, method = "spearman", exact = FALSE)# < 2.2e-16
cor.test(Falphadiv_metadata$Shannon, Falphadiv_metadata$Simpson, method = "spearman", exact = FALSE)# < 2.2e-16

##Multiple test adjustment##
#Make a list of all the p vales from the correlation tests, make it a dataframe, add column of adjusted p values#
alpha_div_p.value=c(7.016e-16,2.2e-16,2.2e-16)
alpha_div_p.value=data.frame(alpha_div_p.value)
alpha_div_p.value$padj <- p.adjust(alpha_div_p.value$alpha_div_p.value, method = "BH")
alpha_div_p.value

##Basic sumamry of Shannon statistics##
mean(Falphadiv_metadata$Shannon)#3.082956
sd(Falphadiv_metadata$Shannon) #0.3439759

ddply(Falphadiv_metadata, .(Diet), summarise, mean=mean(Shannon), sd=sd(Shannon), median=median(Shannon), IQR=IQR(Shannon))

# Diet     mean        sd   median       IQR
# 1   D1 3.235856 0.4680369 3.350047 0.8498954
# 2   D2 2.912338 0.2165244 2.949507 0.2810519
# 3   D3 3.132363 0.3592040 3.129123 0.3132330
# 4   D4 3.051267 0.2226174 3.030987 0.3409434

#Make a simple boxplot to get a visual summary of these statistics#
ggplot(Falphadiv_metadata, aes(x=Diet, y=Shannon, color=Diet)) + 
  geom_boxplot() +
  theme(axis.text.x  = element_text(angle=45, hjust=1))

###TESTS FOR TOTAL AND PAIRWISE SIGNIFICANCE###

#Anova
res.aov <- aov(Shannon ~ Diet, data = Falphadiv_metadata)
summary(res.aov)
#             Df Sum Sq Mean Sq F value  Pr(>F)
# Diet         3  0.671  0.2237   2.013  0.126
# Residuals   44  4.890  0.1111      

#Tukey Test
TukeyHSD(res.aov)
# $Diet
# diff        lwr       upr     p adj
# D2-D1 -0.32351807 -0.6868937 0.0398576 0.0965462
# D3-D1 -0.10349303 -0.4668687 0.2598826 0.8716922
# D4-D1 -0.18458848 -0.5479641 0.1787872 0.5329158
# D3-D2  0.22002503 -0.1433506 0.5834007 0.3799873
# D4-D2  0.13892959 -0.2244461 0.5023053 0.7382166
# D4-D3 -0.08109545 -0.4444711 0.2822802 0.9327881

##Kruskal-Wallis
kruskal.test(Shannon ~ Diet, data = Falphadiv_metadata)
# Kruskal-Wallis chi-squared = 4.9158, df = 3, p-value = 0.1781
##Dunn:
dunnTest(Shannon ~ Diet, data = Falphadiv_metadata, method="bh") 
# Comparison          Z    P.unadj     P.adj
# 1    D1 - D2  2.1141429 0.03450306 0.2070184
# 2    D1 - D3  0.5248907 0.59965920 0.5996592
# 3    D2 - D3 -1.5892523 0.11200345 0.3360103
# 4    D1 - D4  1.0935222 0.27416458 0.5483292
# 5    D2 - D4 -1.0206207 0.30743417 0.4611512
# 6    D3 - D4  0.5686315 0.56960621 0.6835275


###COMBINING STATISTICAL AND ALPHA DIVERSITY TESTING###

##SUPPLEMENTARY FIGURE 1

##Making a more informative boxplot##

#Most of this you've seen before with bar plot#
Shannon_Diet_BP <- ggplot(Falphadiv_metadata, aes(x=Diet, y=Shannon)) + 
  geom_boxplot() +
  xlab("Clam Meal Percentage") +
  annotate("text", x = c(1:4) , y = c(3.85, 3.3, 3.9, 3.5), label = c("a", "a", "a", "a"), size=5)+
  scale_x_discrete("Clam Meal Percentage", labels=c("0%", "10%", "20%", "30%"))+
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=11), legend.key.size = unit(0.7,"cm"), legend.title=element_text(size=13))
Shannon_Diet_BP

tiff('Shannon Diversity by Diet.tiff', units="in", width=10, height=6, res=300)
Shannon_Diet_BP
dev.off()

###BETA DIVERSITY ANALYSIS OTHER TRANSFORMATIONS VERSION###

#Test three different transformations to see which one works best in terms of nomrality and reducing differences between sequencing depths

##Transform data using square root
SRFdata <- transform_sample_counts(Fdata, function(x){x^(1/2)})

##Transform data using DESeq2
#Alternative:
#https://github.com/joey711/phyloseq/issues/299
#https://github.com/joey711/phyloseq/issues/229

#Load alternative gm_mean function to handle zeros based upon github issue shown above#
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# Use `~1` as the experimental design so that the actual design doesn't influence your tranformation.
dds = phyloseq_to_deseq2(Fdata, ~1)

#Calculate geometric means prior to estimate size factors#
geoMeans = apply(counts(dds), 1, gm_mean)
dds = estimateSizeFactors(dds, geoMeans = geoMeans)

#Conduct DESEQ2 test#
dds = DESeq(dds, fitType="local")

#Make a copy of Fdata so you can have a designated DESeq transformed copy
Fdata
DSFdata = Fdata
DSFdata

#Switch the asv table with the DESeq2 transformed data
otu_table(DSFdata) <- otu_table(getVarianceStabilizedData(dds), taxa_are_rows = TRUE)

#Check to see if your basic phyloseq information was not changed
DSFdata

#https://bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis
#https://github.com/joey711/phyloseq/issues/445

#Make any negative values equal to 0 since the more negative they are the more likely they were to be zero over very small and unlikely to affect your results

#Make a copy of DSFdata to manipulate
ZDSFdata <- DSFdata
ZDSFdata

#extract out asv table
DESeq2_otu_table <- as.data.frame(otu_table(ZDSFdata))

#Change all negatives to zero
DESeq2_otu_table[DESeq2_otu_table < 0.0] <- 0.0

#Switch out the asv table in phyloseq object
otu_table(ZDSFdata) <-otu_table(DESeq2_otu_table, taxa_are_rows = TRUE)

#Check to make sure basic phyloseq info did not change
ZDSFdata

#Show how amount of positive numbers changed throughout the transformation 
z <- otu_table(Fdata)
table(as.vector(z) > 0) / prod(dim(z))

# FALSE       TRUE 
# 0.94977227 0.05022773 

z <- otu_table(DSFdata)
table(as.vector(z) > 0) / prod(dim(z))

# FALSE       TRUE 
# 0.94977227 0.05022773 

z <- otu_table(ZDSFdata)
table(as.vector(z) > 0) / prod(dim(z))

# FALSE       TRUE 
# 0.94977227 0.05022773

##Standardize with CSS:
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

#Convert file types to use in metagenomeSeq
MGS <- phyloseq_to_metagenomeSeq(Fdata) 

#Perform normalization following: https://bioconductor.org/packages/release/bioc/vignettes/metagenomeSeq/inst/doc/metagenomeSeq.pdf
CSSp <- cumNormStatFast(MGS)
CSSp
MGS <- cumNorm(MGS, p = CSSp)

#Store post-transformation abundance table
normmybiom <- MRcounts(MGS, norm = T)

#Create copy of filtered data to be modified
MGSdata = Fdata
MGSdata

#Switch out old ASV table for transformed one
otu_table(MGSdata) <-otu_table(normmybiom, taxa_are_rows = TRUE)
MGSdata

#Check your library size differences by dividing the sequences sum of the sample with the most sequences by the one with the least sequences typically want it to be <10x difference


max(sample_sums(Fdata))/min(sample_sums(Fdata))#9.257675
max(sample_sums(SRFdata))/min(sample_sums(SRFdata)) #5.699979
max(sample_sums(ZDSFdata))/min(sample_sums(ZDSFdata)) #3.864704
max(sample_sums(MGSdata))/min(sample_sums(MGSdata)) #2.815309

#Check to see if they are statistically normal, if >0.05 then it is likely normal
shapiro.test(sample_sums(Fdata)) #0.3178
shapiro.test(sample_sums(SRFdata)) #0.2794
shapiro.test(sample_sums(ZDSFdata)) #0.02869
shapiro.test(sample_sums(MGSdata)) #0.01584

#Visual representation of the sample sums
hist(sample_sums(Fdata))
hist(sample_sums(SRFdata))
hist(sample_sums(ZDSFdata))
hist(sample_sums(MGSdata))

#Typically want to ratio of high/low to be less than 10 and a more normal distibution, thus the square root transfomration was chosen

###BETA DIVERSITY ANALYSIS SQUARE ROOT VERSION###

##Export ESV table for PRIMER7##

# Extract abundance matrix from the phyloseq object#
QIIME2_WS_Bio = as(otu_table(SRFdata), "matrix")

# Coerce to data.frame#
QIIME2_WS_Bio = as.data.frame(QIIME2_WS_Bio)

#Make a txt file#
write.table(QIIME2_WS_Bio, file = 'QIIME2_WS_Bio.txt', sep = "\t")

##PCoA Analysis##

#Make a dissimilarity/similarity matrix, by default Bray Curtis is used#
SRFdataBrayPCoAO <- ordinate(SRFdata,"PCoA")

#With Tanks
plot_ordination(SRFdata, SRFdataBrayPCoAO, color="Tanks", shape="Diet")+
  geom_point(size=7, alpha=0.75)+
  ggtitle("PCOA with Bray-Curtis")+
  scale_colour_manual(name="Tank", values=c("blue", "orange", "green", "purple", "lightblue", "lightsalmon", "forestgreen", "orchid", "greenyellow", "navyblue", "dodgerblue", "sienna", "olivedrab", "chocolate", "darkorchid", "thistle"))+
  scale_shape_manual(name="Diet", values = c(15,16,17,18))

#Without Tanks

##SUPPLEMENTARY FIGURE 2

Diet_PCoA <- plot_ordination(SRFdata, SRFdataBrayPCoAO, color="Diet")+
  geom_point(size=7, alpha=0.75)+
  #ggtitle("PCOA with Bray-Curtis")+
  scale_colour_manual(name="Clam Meal Percentage", labels = c("0%", "10%", "20%", "30%"), values=c("white", "#70AD46", "#ED7D31", "#5B9BD5"))+
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=11), legend.key.size = unit(0.7,"cm"), legend.title=element_text(size=13))
Diet_PCoA

#remove initial layer
Diet_PCoA$layers
Diet_PCoA$layers[[1]] <- NULL
Diet_PCoA

tiff('PCoA by Diet.tiff', units="in", width=10, height=6, res=300)
Diet_PCoA
dev.off()

###VENN DIAGRMS###

#Subset samples associated with control and remove any ASVs that are no longer represented
D1_Fdata <- subset_samples(Fdata, Diet=="D1")
D1_Fdata #988 taxa
D1_Fdata <- filter_taxa(D1_Fdata, function(x) sum(x) >0, TRUE)
D1_Fdata #480 taxa

#Extract out the asv table and remove the taxonomic labels
D1_Fdata_taxa <- as.data.frame(tax_table(D1_Fdata))

D1_Fdata_taxa <- data.frame(lapply(D1_Fdata_taxa, function(x) {gsub("d__", "", x)}))
D1_Fdata_taxa <- data.frame(lapply(D1_Fdata_taxa, function(x) {gsub("p__", "", x)}))
D1_Fdata_taxa <- data.frame(lapply(D1_Fdata_taxa, function(x) {gsub("c__", "", x)}))
D1_Fdata_taxa <- data.frame(lapply(D1_Fdata_taxa, function(x) {gsub("o__", "", x)}))
D1_Fdata_taxa <- data.frame(lapply(D1_Fdata_taxa, function(x) {gsub("f__", "", x)}))
D1_Fdata_taxa <- data.frame(lapply(D1_Fdata_taxa, function(x) {gsub("g__", "", x)}))
D1_Fdata_taxa <- data.frame(lapply(D1_Fdata_taxa, function(x) {gsub("s__", "", x)}))

#Create a table agglomerated to genus level
D1_Fdata_taxa_genus <- D1_Fdata_taxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, .keep_all=TRUE)

#Fuse the entire taxonomy together into one column, cannot just compare the genus level because multiple occurrences of "same" genus name, like g__uncultured, with different preceding taxonomy strings
D1_Fdata_taxa_genus$Taxonomy <- paste(D1_Fdata_taxa_genus$Kingdom,D1_Fdata_taxa_genus$Phylum,D1_Fdata_taxa_genus$Class,D1_Fdata_taxa_genus$Order,D1_Fdata_taxa_genus$Family,D1_Fdata_taxa_genus$Genus,sep="-")

#Extract out the full taxonomy list column as a list
D1_Fdata_genera_list <- as.list(D1_Fdata_taxa_genus$Taxonomy)#113


#Do the same for the 30% clam meal samples
D4_Fdata <- subset_samples(Fdata, Diet=="D4")
D4_Fdata #988 taxa
D4_Fdata <- filter_taxa(D4_Fdata, function(x) sum(x) >0, TRUE)
D4_Fdata #371 taxa

#Extract out the asv table and remove the taxonomic labels
D4_Fdata_taxa <- as.data.frame(tax_table(D4_Fdata))

D4_Fdata_taxa <- data.frame(lapply(D4_Fdata_taxa, function(x) {gsub("d__", "", x)}))
D4_Fdata_taxa <- data.frame(lapply(D4_Fdata_taxa, function(x) {gsub("p__", "", x)}))
D4_Fdata_taxa <- data.frame(lapply(D4_Fdata_taxa, function(x) {gsub("c__", "", x)}))
D4_Fdata_taxa <- data.frame(lapply(D4_Fdata_taxa, function(x) {gsub("o__", "", x)}))
D4_Fdata_taxa <- data.frame(lapply(D4_Fdata_taxa, function(x) {gsub("f__", "", x)}))
D4_Fdata_taxa <- data.frame(lapply(D4_Fdata_taxa, function(x) {gsub("g__", "", x)}))
D4_Fdata_taxa <- data.frame(lapply(D4_Fdata_taxa, function(x) {gsub("s__", "", x)}))

#Create a table agglomerated to genus level
D4_Fdata_taxa_genus <- D4_Fdata_taxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, .keep_all=TRUE)

#Fuse the entire taxonomy together into one column, cannot just compare the genus level because multiple occurrences of "same" genus name, like g__uncultured, with different preceding taxonomy strings
D4_Fdata_taxa_genus$Taxonomy <- paste(D4_Fdata_taxa_genus$Kingdom,D4_Fdata_taxa_genus$Phylum,D4_Fdata_taxa_genus$Class,D4_Fdata_taxa_genus$Order,D4_Fdata_taxa_genus$Family,D4_Fdata_taxa_genus$Genus,sep="-")

#Extract out the full taxonomy list column as a list
D4_Fdata_genera_list <- as.list(D4_Fdata_taxa_genus$Taxonomy)#108

#Combine the list into a list of lists
D1vsD4_Genera_Lists <- list('0% Clam Meal (113)' = D1_Fdata_genera_list,
                            '30% Clam Meal (108)' = D4_Fdata_genera_list)

#Create a Venn Diagram comparing Genera and save it
ggvenn(D1vsD4_Genera_Lists, c("0% Clam Meal (113)", "30% Clam Meal (108)"), fill_color = c("white", "#5B9BD5"), show_percentage = FALSE)+
  ggtitle("0% vs 30% Clam Feed Genera") +
  theme(plot.title = element_text(hjust = 0.5))

tiff('D1 vs D4 Clam Feed Genera venn diagram.tiff', units="in", width=10, height=6, res=300)
#plot
dev.off()


###TAXONOMIC BAR PLOTS###
# get abundance in relative percentage#
phy <- transform_sample_counts(Fdata, function(x) 100*x/sum(x))

##Create table ready for making stacked bar graph for Phyla >1% by Diet##

# agglomerate taxa such that the ASVs with the same upper taxonomic ranks are merged at whatever rank you choose#
glom <- tax_glom(phy, taxrank = 'Phylum', NArm=FALSE)

# Create dataframe from phyloseq object, just melts the phyloseq object to a full table#
dat <- psmelt(glom)

# convert Phylum to a character vector from a factor because R#
dat$Phylum <- as.character(dat$Phylum)

# group dataframe by Phylum, calculate mean rel. abundance#
means <- ddply(dat, ~Phylum, function(x) c(mean=mean(x$Abundance)))

# find Phyla whose rel. abund. is less than 1%#
Other <- means[means$mean <= 1,]$Phylum

# change their name to "Other Prokaryotes"#
dat[dat$Phylum %in% Other,]$Phylum <- 'Other Prokaryotes'

#remove all Phylums labeled Other Prokaryotes#
dat <-dat[!dat$Phylum=="Other Prokaryotes",]

#remove unnessary columns#
dat <- subset(dat, select=c(Diet, Abundance, Phylum))

#Summarize based upon target parameter#
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Diet","Phylum"), na.rm = TRUE)

#Arrange by Abundance#
dat <- arrange(dat, Abundance)

#Create a table that is the leftover Prokaryotes#
Abundance <- ddply(dat, ~Diet, function(x) c(Abundance=100-sum(x$Abundance)))

#Add a column labeling the leftover Prokaryotes#
Abundance$Phylum<- "Other Prokaryotes"

#remove unnessary columns#
dat <- subset(dat, select=c(Diet, Abundance, Phylum))

#combine with original table#
Phylum_Fdata_Diet <- rbind(dat, Abundance)

##Create table ready for making stacked bar graph for Classes <0.5%##
# agglomerate taxa
glom <- tax_glom(phy, taxrank = 'Class')
# create dataframe from phyloseq object
dat <- psmelt(glom)
# convert Class to a character vector from a factor because R
dat$Class <- as.character(dat$Class)
# group dataframe by Class, calculate mean rel. abundance
means <- ddply(dat, ~Class, function(x) c(mean=mean(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
Other <- means[means$mean <= 1,]$Class
# change their name to "Other Prokaryotes"
dat[dat$Class %in% Other,]$Class <- 'Other Prokaryotes'
#remove all Genera labeled Other Prokaryotes
dat <-dat[!dat$Class=="Other Prokaryotes",]
#remove unnessary columns
dat <- subset(dat, select=c(Diet, Abundance, Class))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Diet","Class"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#Create a table that is the leftover Prokaryotes
Abundance <- ddply(dat, ~Diet, function(x) c(Abundance=100-sum(x$Abundance)))
#Add a column labeling the leftover Prokaryotes
Abundance$Class<- "Other Prokaryotes"
#remove unnessary columns
dat <- subset(dat, select=c(Diet, Abundance, Class))
#combine with original table
Class_Fdata_Diet <- rbind(dat, Abundance)

##Create table ready for making stacked bar graph for Orders <0.5%##
# agglomerate taxa
glom <- tax_glom(phy, taxrank = 'Order')
# create dataframe from phyloseq object
dat <- psmelt(glom)
# convert Order to a character vector from a factor because R
dat$Order <- as.character(dat$Order)
# group dataframe by Order, calculate mean rel. abundance
means <- ddply(dat, ~Order, function(x) c(mean=mean(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
Other <- means[means$mean <= 1,]$Order
# change their name to "Other Prokaryotes"
dat[dat$Order %in% Other,]$Order <- 'Other Prokaryotes'
#remove all Genera labeled Other Prokaryotes
dat <-dat[!dat$Order=="Other Prokaryotes",]
#remove unnessary columns
dat <- subset(dat, select=c(Diet, Abundance, Order))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Diet","Order"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#Create a table that is the leftover Prokaryotes
Abundance <- ddply(dat, ~Diet, function(x) c(Abundance=100-sum(x$Abundance)))
#Add a column labeling the leftover Prokaryotes
Abundance$Order<- "Other Prokaryotes"
#remove unnessary columns
dat <- subset(dat, select=c(Diet, Abundance, Order))
#combine with original table
Order_Fdata_Diet <- rbind(dat, Abundance)

##Create table ready for making stacked bar graph for Families <0.5%##
# agglomerate taxa
glom <- tax_glom(phy, taxrank = 'Family')
# create dataframe from phyloseq object
dat <- psmelt(glom)
# convert Family to a character vector from a factor because R
dat$Family <- as.character(dat$Family)
# group dataframe by Family, calculate mean rel. abundance
means <- ddply(dat, ~Family, function(x) c(mean=mean(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
Other <- means[means$mean <= 1,]$Family
# change their name to "Other Prokaryotes"
dat[dat$Family %in% Other,]$Family <- 'Other Prokaryotes'
#remove all Genera labeled Other Prokaryotes
dat <-dat[!dat$Family=="Other Prokaryotes",]
#remove unnessary columns
dat <- subset(dat, select=c(Diet, Abundance, Family))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Diet","Family"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#Create a table that is the leftover Prokaryotes
Abundance <- ddply(dat, ~Diet, function(x) c(Abundance=100-sum(x$Abundance)))
#Add a column labeling the leftover Prokaryotes
Abundance$Family<- "Other Prokaryotes"
#remove unnessary columns
dat <- subset(dat, select=c(Diet, Abundance, Family))
#combine with original table
Family_Fdata_Diet <- rbind(dat, Abundance)

##Create table ready for making stacked bar graph for Genera <0.5%##
# agglomerate taxa
glom <- tax_glom(phy, taxrank = 'Genus')
# create dataframe from phyloseq object
dat <- psmelt(glom)
#Combine Family and Genus columns#
dat<- unite(dat, Genus, Family:Genus, sep=';')
# convert Genus to a character vector from a factor because R
dat$Genus <- as.character(dat$Genus)
# group dataframe by Genus, calculate mean rel. abundance
means <- ddply(dat, ~Genus, function(x) c(mean=mean(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
Other <- means[means$mean <= 1,]$Genus
# change their name to "Other Prokaryotes"
dat[dat$Genus %in% Other,]$Genus <- 'Other Prokaryotes'
#remove all Genera labeled Other Prokaryotes
dat <-dat[!dat$Genus=="Other Prokaryotes",]
#remove unnessary columns
dat <- subset(dat, select=c(Diet, Abundance, Genus))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Diet","Genus"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#Create a table that is the leftover Prokaryotes
Abundance <- ddply(dat, ~Diet, function(x) c(Abundance=100-sum(x$Abundance)))
#Add a column labeling the leftover Prokaryotes
Abundance$Genus<- "Other Prokaryotes"
#remove unnessary columns
dat <- subset(dat, select=c(Diet, Abundance, Genus))
#combine with original table
Genus_Fdata_Diet <- rbind(dat, Abundance)


##Make a series of graphs using ggplot based on above tables##

#Make Diet Phylum graph#

spatial_plot_Phylum_Fdata_Diet <- ggplot(data=Phylum_Fdata_Diet, aes(x=Diet, y=Abundance, fill=Phylum)) + 
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("Phylum", values=c("deeppink", "midnightblue", "blue", "rosybrown", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "lightblue", "turquoise", "firebrick3"))+
  xlab("Site by Survey") +
  ylab("Percentage") +
  ggtitle("Phyla by Site by Survey")+
  theme(plot.title = element_text(hjust = 0.5)) 
spatial_plot_Phylum_Fdata_Diet

##Save Graph for Publication##
tiff('Phylum by Site by Survey.tiff', units="in", width=10, height=6, res=300)
spatial_plot_Phylum_Fdata_Diet
dev.off()

##Make Class Graph##

spatial_plot_Class_Fdata_Diet <- ggplot(data=Class_Fdata_Diet, aes(x=Diet, y=Abundance, fill=Class)) + 
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("Class", values=c("deeppink", "midnightblue", "blue", "rosybrown", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "lightblue", "turquoise", "firebrick3"))+
  xlab("Site by Survey") +
  ylab("Percentage") +
  ggtitle("Phyla by Site by Survey")+
  theme(plot.title = element_text(hjust = 0.5)) 
spatial_plot_Class_Fdata_Diet

##Make Order Graph##

spatial_plot_Order_Fdata_Diet <- ggplot(data=Order_Fdata_Diet, aes(x=Diet, y=Abundance, fill=Order)) + 
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("Order", values=c("deeppink", "midnightblue", "blue", "rosybrown", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "lightblue", "turquoise", "firebrick3"))+
  xlab("Site by Survey") +
  ylab("Percentage") +
  ggtitle("Phyla by Site by Survey")+
  theme(plot.title = element_text(hjust = 0.5)) 
spatial_plot_Order_Fdata_Diet

##Make Family Graph##

spatial_plot_Family_Fdata_Diet <- ggplot(data=Family_Fdata_Diet, aes(x=Diet, y=Abundance, fill=Family)) + 
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("Family", values=c("deeppink", "midnightblue", "blue", "rosybrown", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "lightblue", "turquoise", "firebrick3"))+
  xlab("Site by Survey") +
  ylab("Percentage") +
  ggtitle("Phyla by Site by Survey")+
  theme(plot.title = element_text(hjust = 0.5)) 
spatial_plot_Family_Fdata_Diet

##Make Genus Graph##

##SUPPLEMENTARY FIGURE 3

spatial_plot_Genus_Fdata_Diet <- ggplot(data=Genus_Fdata_Diet, aes(x=Diet, y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  #theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("Family;Genus", values=c("deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "lightblue", "orangered", "firebrick3", "darkseagreen3", "navy", "black", "seagreen", "turquoise", "rosybrown", "white", "darkviolet", "wheat", "plum3"))+
  #ggtitle("Genera by Diet") +
  xlab("Clam Meal Percentage") +
  scale_x_discrete(labels=c("0%", "10%", "20%", "30%"))+
  ylab("Family;Genus Percentage") +
  #theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=11), legend.key.size = unit(0.7,"cm"), legend.title=element_text(size=13))
spatial_plot_Genus_Fdata_Diet

tiff('Genera by Diet.tiff', units="in", width=10, height=6, res=300)
spatial_plot_Genus_Fdata_Diet
dev.off()

###CREATE TABLES OF TAXONOMIC LEVELS TO GET SUMMARY STATISTICS FOR EACH ENTRY

##Create table with all Phyla##
# agglomerate taxa
glom <- tax_glom(phy, taxrank = 'Phylum')
# create dataframe from phyloseq object
dat <- psmelt(glom)
# convert Phylum to a character vector from a factor because R
dat$Phylum <- as.character(dat$Phylum)
#remove unnessary columns
dat <- subset(dat, select=c(Abundance, Phylum))
#Summarize based upon target parameter
Phylum_Summary <- summarySE(data=dat, measurevar="Abundance", groupvars="Phylum", na.rm = TRUE)

##Create table with all Genera##
# agglomerate taxa
glom <- tax_glom(phy, taxrank = 'Genus')
# create dataframe from phyloseq object
dat <- psmelt(glom)
#Combine Family and Genus columns#
dat<- unite(dat, Genus, Family:Genus, sep=';')
# convert Genus to a character vector from a factor because R
dat$Genus <- as.character(dat$Genus)
#remove unnessary columns
dat <- subset(dat, select=c(Abundance, Genus))
#Summarize based upon target parameter
Genus_Summary <- summarySE(data=dat, measurevar="Abundance", groupvars="Genus", na.rm = TRUE)

###CREATE TABLE FOR ORGANISM BY ORGANISM COMPARISONS BETWEEN SAMPLE TYPES

#Phylum Level

##Create table ready for making stacked bar graph for Phylumes <0.5%##
# agglomerate taxa
glom <- tax_glom(phy, taxrank = 'Phylum')
# create dataframe from phyloseq object
dat <- psmelt(glom)
# convert Phylum to a character vector from a factor because R
dat$Phylum <- as.character(dat$Phylum)
#remove unnessary columns
dat <- subset(dat, select=c(Diet, Abundance, Phylum))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Diet","Phylum"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#remove unnessary columns
dat <- subset(dat, select=c(Diet, Abundance, Phylum))
#combine with original table
All_Phylum_Fdata_Diet <- dat

#Switch the structure of dataframe from long to wide
wide_All_Phylum_Fdata_Diet <- spread(All_Phylum_Fdata_Diet, Diet, Abundance) 

#Add a new column that is the Diet 1 abundance of a particular genera divided by the Diet 4 abundance of that same genera
wide_All_Phylum_Fdata_Diet$D1D4 <- wide_All_Phylum_Fdata_Diet$D1 / wide_All_Phylum_Fdata_Diet$D4 

#Same thing except that its D4 divided by D1
wide_All_Phylum_Fdata_Diet$D4D1 <- wide_All_Phylum_Fdata_Diet$D4 / wide_All_Phylum_Fdata_Diet$D1

#Add a new column that is the D1 abundance of a particular genera divided by the D3 abundance of that same genera
wide_All_Phylum_Fdata_Diet$D1D3 <- wide_All_Phylum_Fdata_Diet$D1 / wide_All_Phylum_Fdata_Diet$D3 

#Same thing except that its D3 divided by D1
wide_All_Phylum_Fdata_Diet$D3D1 <- wide_All_Phylum_Fdata_Diet$D3 / wide_All_Phylum_Fdata_Diet$D1

#Add a new column that is the D2 abundance of a particular genera divided by the D3 abundance of that same genera
wide_All_Phylum_Fdata_Diet$D2D3 <- wide_All_Phylum_Fdata_Diet$D2 / wide_All_Phylum_Fdata_Diet$D3 

#Same thing except that its D3 divided by D2
wide_All_Phylum_Fdata_Diet$D3D2 <- wide_All_Phylum_Fdata_Diet$D3 / wide_All_Phylum_Fdata_Diet$D2


##Genus Level

##Create table ready for making stacked bar graph for Genuses <0.5%##
# agglomerate taxa
glom <- tax_glom(phy, taxrank = 'Genus')
# create dataframe from phyloseq object
dat <- psmelt(glom)
#Combine Family and Genus columns#
dat<- unite(dat, Genus, Family:Genus, sep=';')
# convert Genus to a character vector from a factor because R
dat$Genus <- as.character(dat$Genus)
#remove unnessary columns
dat <- subset(dat, select=c(Diet, Abundance, Genus))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Diet","Genus"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#remove unnessary columns
dat <- subset(dat, select=c(Diet, Abundance, Genus))
#combine with original table
All_Genus_Fdata_Diet <- dat


#Switch the structure of dataframe from long to wide
wide_All_Genus_Fdata_Diet <- spread(All_Genus_Fdata_Diet, Diet, Abundance) 

#Add a new column that is the D1 abundance of a particular genera divided by the D4 abundance of that same genera
wide_All_Genus_Fdata_Diet$D1D4 <- wide_All_Genus_Fdata_Diet$D1 / wide_All_Genus_Fdata_Diet$D4 

#Same thing except that its D4 divided by D1
wide_All_Genus_Fdata_Diet$D4D1 <- wide_All_Genus_Fdata_Diet$D4 / wide_All_Genus_Fdata_Diet$D1

#Add a new column that is the D1 abundance of a particular genera divided by the D3 abundance of that same genera
wide_All_Genus_Fdata_Diet$D1D3 <- wide_All_Genus_Fdata_Diet$D1 / wide_All_Genus_Fdata_Diet$D3 

#Same thing except that its D3 divided by D1
wide_All_Genus_Fdata_Diet$D3D1 <- wide_All_Genus_Fdata_Diet$D3 / wide_All_Genus_Fdata_Diet$D1

#Add a new column that is the D2 abundance of a particular genera divided by the D3 abundance of that same genera
wide_All_Genus_Fdata_Diet$D2D3 <- wide_All_Genus_Fdata_Diet$D2 / wide_All_Genus_Fdata_Diet$D3 

#Same thing except that its D3 divided by D2
wide_All_Genus_Fdata_Diet$D3D2 <- wide_All_Genus_Fdata_Diet$D3 / wide_All_Genus_Fdata_Diet$D2

###DIFFERENTIAL ABUNDANCE

##ASV Level

##Subset out the 0% and 30% clam meal samples

#uf means it that the subseted phyloseq has not been filtered of ASVs that now have 0 sequences, this is done so that these phyloseqs can be combined, otherwise their trees would not match if they had been filtered again

D4_Fdata_uf <- subset_samples(Fdata, Diet=="D4")
D1_Fdata_uf <- subset_samples(Fdata, Diet=="D1")
D14_Fdata =merge_phyloseq(D4_Fdata_uf,D1_Fdata_uf)

# Convert the phyloseq object to a DESeq2 one, this time using Diet as the variable instead of 1 as used above for the transformation test
dds = phyloseq_to_deseq2(D14_Fdata, ~Diet)

#Calculate geometric means prior to estimate size factors#
geoMeans = apply(counts(dds), 1, gm_mean)
dds = estimateSizeFactors(dds, geoMeans = geoMeans)

#Conduct DESEQ2 test#
dds = DESeq(dds, fitType="local")

#Explore the results#
res_asv = results(dds)
res_asv = res_asv[order(res_asv$padj, na.last=NA), ]
alpha = 0.05
sig_res_asv = res_asv[(res_asv$padj < alpha), ]

#Make dataframe with taxanomy added in#
sig_res_taxo_asv = cbind(as(sig_res_asv, "data.frame"), as(tax_table(D14_Fdata)[rownames(sig_res_asv), ], "matrix"))

#Make dataframe with asvs added in from all sites#
sig_res_taxo_seqs_asv = cbind(as(sig_res_taxo_asv, "data.frame"), as(otu_table(D14_Fdata)[rownames(sig_res_taxo_asv), ], "matrix"))

#Make rownames an actual column and remove old rownames#
sig_res_taxo_seqs_asv <- cbind(ESV.ID = rownames(sig_res_taxo_seqs_asv), sig_res_taxo_seqs_asv)
rownames(sig_res_taxo_seqs_asv) <-NULL

#Make a txt file of results for all sites#
write.csv(as.data.frame(sig_res_taxo_seqs_asv), file="sig_res_taxo_seqs_asv.csv")

#Detemine which subcategory is negative or positive#
res_asv
#res_asv: D4 (+) vs D1 (-)#

#Determine number of Indicators with padj >0.05 in each subcagetory#
#Number of D4 indicators# 
length(which(sig_res_taxo_seqs_asv$log2FoldChange > 0)) #0
#Number of D1 indicators#
length(which(sig_res_taxo_seqs_asv$log2FoldChange < 0)) #2

##Genus Level

##Upload Genus level files from QIIME2
genus_aa_table = read.csv(file.choose(), stringsAsFactors = TRUE, row.names=1)
colnames(genus_aa_table) <- gsub("^X", "",  colnames(genus_aa_table))
genus_OTU = otu_table(genus_aa_table, taxa_are_rows = TRUE)
sample_names(genus_OTU)

sample_names(META)

genus_taxa_table = read.csv(file.choose(), row.names = 1)
genus_taxa_table = as.matrix(genus_taxa_table)
genus_TAX = tax_table(genus_taxa_table)

genus <- merge_phyloseq(genus_OTU, genus_TAX, META)
genus #272

Fgenus = filter_taxa(genus, function(x) sum(x) >20, TRUE)
Fgenus #177

##Subset out the 0% and 30% clam meal samples
D4_Fgenus_uf <- subset_samples(Fgenus, Diet=="D4")
D1_Fgenus_uf <- subset_samples(Fgenus, Diet=="D1")
D14_Fgenus =merge_phyloseq(D4_Fgenus_uf,D1_Fgenus_uf)


# Convert the phyloseq object to a DESeq2 one
dds = phyloseq_to_deseq2(D14_Fgenus, ~Diet)

#Calculate geometric means prior to estimate size factors#
geoMeans = apply(counts(dds), 1, gm_mean)
dds = estimateSizeFactors(dds, geoMeans = geoMeans)

#Conduct DESEQ2 test#
dds = DESeq(dds, fitType="local")

#Explore the results#
res_genus = results(dds)
res_genus = res_genus[order(res_genus$padj, na.last=NA), ]
alpha = 0.05
sig_res_genus = res_genus[(res_genus$padj < alpha), ]

#Make dataframe with taxonomy added in#
sig_res_taxo_genus = cbind(as(sig_res_genus, "data.frame"), as(tax_table(D14_Fgenus)[rownames(sig_res_genus), ], "matrix"))

#Fails because there are no diferentially abundant genera