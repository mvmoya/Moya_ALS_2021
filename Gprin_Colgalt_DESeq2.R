##########################################################################################################
###                                                                                                    ###
### This is the script I used to load and analyze the files from the May-June 2019 Gprin3 and Colgalt2 ###
### IP experiments. It has been optimized to work for DESeq2 input and output.                         ###
###                                                                                                    ###
##########################################################################################################

# Here, we load the libraries that we'll need to run all of the analysis. Some of these can be deprecated

library(DESeq2)
library(edgeR)
library(reshape2)
library(ggplot2)
library(Hmisc)
library(Heatplus)
library(RColorBrewer)
library(pheatmap)
library(ggfortify)
library(dplyr)
library(stringr)
library(ggrepel)
library(topGO)
library(clipr)

# Here, set the working directory to wherever the counts files are all housed:

setwd('My/working/directory')

# This line gets a list of the files in the directory if they have the word "ID" in them:

sampleFiles <- grep("(Gprin3|Colgalt2|Glt25d2).*(Counts)", list.files(getwd()), value = TRUE)

# Here, we create a list of what we want the sample names to be. As written, this has to match the
# order of the files on the server. Sorry for the brute force approach:

sampleCondition <- c('Colgalt2_ALS_WT_Pre_IP_1', 'Colgalt2_ALS_WT_Pre_IP_2', 'Colgalt2_ALS_WT_Pre_IP_3',
                     'Colgalt2_ALS_SOD_Pre_IP_1', 'Colgalt2_ALS_SOD_Pre_IP_2',
                     'Colgalt2_ALS_SOD_Pre_IP_3', 'Colgalt2_ALS_WT_Post_IP_1', 
                     'Colgalt2_ALS_WT_Post_IP_3',
                     'Colgalt2_ALS_SOD_Post_IP_1', 'Colgalt2_ALS_SOD_Post_IP_2',
                     'Colgalt2_ALS_SOD_Post_IP_3', 'Colgalt2_ALS_WT_Post_IP_4',
                     'Colgalt2_ALS_SOD_Post_IP_4', 'Colgalt2_ALS_WT_Pre_Input',
                     'Colgalt2_ALS_SOD_Pre_Input', 'Colgalt2_ALS_WT_Post_Input',
                     'Gprin3_Norm_M1_IP_1', 'Gprin3_Norm_M1_IP_2', 'Gprin3_Norm_M1_IP_3',
                     'Gprin3_Norm_M1_Input_1', 'Gprin3_Norm_M1_Input_2', 'Gprin3_Norm_M1_Input_3',
                     'Gprin3_Norm_Striat_IP_1', 'Gprin3_Norm_Striat_IP_2', 'Gprin3_Norm_Striat_IP_3',
                     'Colgalt2_Norm_M1_IP_1', 'Colgalt2_Norm_M1_IP_2', 'Colgalt2_Norm_M1_Input_1',
                     'Colgalt2_Norm_M1_Input_2', 'Colgalt2_Norm_M1_IP_3', 'Colgalt2_Norm_M1_Input_3',
                     'Colgalt2_Norm_Whole_IP_1', 'Colgalt2_Norm_Whole_IP_2',
                     'Gprin3_Norm_Whole_IP_1', 'Gprin3_Norm_Whole_IP_2',
                     'Gprin3_ALS_SOD_Pre_IP_1', 'Gprin3_ALS_SOD_Pre_IP_2',
                     'Gprin3_ALS_SOD_Pre_IP_3', 'Gprin3_ALS_SOD_Post_IP_1',
                     'Gprin3_ALS_SOD_Post_IP_2', 'Gprin3_ALS_SOD_Post_IP_3',
                     'Gprin3_ALS_SOD_Early_IP_1', 'Gprin3_ALS_SOD_Early_IP_2',
                     'Gprin3_ALS_SOD_Early_IP_3', 'Gprin3_ALS_WT_Pre_IP_1', 
                     'Gprin3_ALS_WT_Pre_IP_2', 'Gprin3_ALS_WT_Pre_IP_3', 'Gprin3_ALS_WT_Post_IP_1',
                     'Gprin3_ALS_WT_Post_IP_2', 'Gprin3_ALS_WT_Post_IP_3',
                     'Gprin3_ALS_WT_Early_IP_1', 'Gprin3_ALS_WT_Early_IP_3')

# Now we make a reference table for the sample names from above and the corresponding filename:

sampleTable <- data.frame(sampleName = sampleFiles, fileName = sampleFiles, condition = sampleCondition)

# Here we read the first file in the list and seed a data frame from it. I've changed the header
# names to be easier to reference, and have removed rows where gene names are NA values:

data.num <- read.delim(toString(sampleTable[1, c('fileName')]), header = TRUE)
data.num <- data.num[,2:3] # For HPC output
#data.num <- data.num[1:32658,] # For Heintz pipeline output
data.num <- data.num[complete.cases(data.num),]
colnames(data.num) <- c('genes', sampleCondition[1])

# Now, we add the remaining sample files to this table by merging the tables by gene name.
# This part loops through every filename specified in our reference table above:

for (i in 2:nrow(sampleTable)) {
  new_cols <- read.delim(toString(sampleTable[i,"fileName"]), header = TRUE)
  new_cols <- new_cols[,2:3] # For HPC output
  #new_cols <- new_cols[1:32658,] # For Heintz pipeline output
  colnames(new_cols) <- c('genes', sampleCondition[i])
  data.num <- merge(data.num, new_cols, by = "genes")
}

# Now we make the rownames for the data frame into the gene names:

rownames(data.num) <- data.num[,1]
data.num[,1] <- NULL

# Here, we create the experimental information table. Using regular expressions, I make a list
# that will tell DESeq if our sample was a Gprin3 or Colgalt2 IP or Input sample:

#data.num <- subset(data.num, select=-c(Colgalt2_ALS_WT_Post_IP_2, Colgalt2_ALS_SOD_Post_Input)) #Colgalt2_nIP_2, Colgalt2_nIP_3, 
#sample_type <- as.vector(sub("(.*_.*)_.*", "\\1", colnames(data.num)))

cell <- sub(".*(Gprin3|Colgalt2).*", "\\1", sampleCondition)
type <- sub(".*(ALS|Norm).*", "\\1", sampleCondition)
disease <- sub(".*(WT|SOD|Norm).*", "\\1", sampleCondition)
disease <- sub("Norm", "WT", disease)
ip_input <- sub(".*(IP|Input).*", "\\1", sampleCondition)
age <- sub(".*(Pre|Post|Early|Norm).*", "\\1", sampleCondition)
age <- sub("Norm", "Pre", age)
when <- sub(".*(M1|Striat).*", "\\1", sampleCondition)
when <- sub(".*(ALS|Whole).*", "Whole", when)

coldata <- data.frame(cell, type, disease, ip_input, age, when)
colnames(coldata) <- factor(c('cell', 'type', 'disease', 'ip_input', 'age', 'when'))
# Now, we create the DESeq object, specifying that we want sample_type to be what groups our samples:

ddsPre <- DESeqDataSetFromMatrix(countData = data.num, colData = coldata, design = ~ cell + when)
dds <- DESeq(ddsPre)

# These lines create a table of normalized counts from the DESeq object above:

rld <- rlog(dds, blind = FALSE)
norm_counts <- 2**assay(rld)

# (Optional) And we can print out counts values for specific genes of interest with this command:

genes_of_interest <- c('Crym', 'Nefh', 'Gprin3', 'Colgalt2', 'Fezf2', 'Bcl11b', 'Pcp4', 
                       'Slco2a1', 'Npnt')
glial_genes <- c('Aldh1l1', 'Aif1', 'Gfap', 'Pvalb')
norm_counts[c(genes_of_interest, glial_genes),grep('IP', colnames(norm_counts))]

# Here, we generate a PCA of our different sample types. The color levels should match the groups
# This is a rather convoluted way of doing it, but gives you more control over which PCs are plotted

vsd_of_interest <- rld
color_scheme <- factor(as.numeric(factor(when))*(as.numeric(factor(cell))+1))

rv <- rowVars(assay(vsd_of_interest))
select <- order(rv, decreasing=TRUE)#[seq_len(500)]
pca <- prcomp(t(assay(vsd_of_interest)[select,]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
d <- data.frame(PC1 = pca$x[, 1], PC3 = pca$x[, 3], group = cell, name = colnames(data.num))
ggplot(data = d, aes_string(x = 'PC1', y = 'PC3', color = color_scheme, label = 'name')) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = name), size = 2, box.padding = 0.15, 
                  point.padding = 0.5, segment.color = 'grey50') +
  #geom_label_repel(aes(label = name), size = 2, box.padding = 0.35, 
  #                 point.padding = 0.5, segment.color = 'grey50') +
  xlab(paste('PC1: ', round(percentVar[1] * 100), '% variance')) +
  ylab(paste('PC3: ', round(percentVar[3] * 100), '% variance')) +
  coord_fixed() +
  scale_fill_discrete(name='Group', labels=c('Colgalt2 Old', 'Gprin3 Old', 'Colgalt2 New', 'Gprin3 New'))+
  ggtitle('PCA of IP and Input samples')

View(pca$rotation)

#print comp, rotation, pc loading print

# Here, we perform heirarchical clustering of our samples based on gene expression:

sampleDists <- dist(t(assay(vsd_of_interest)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd_of_interest)
colnames(sampleDistMatrix) <- colnames(vsd_of_interest)
pheatmap(sampleDistMatrix, clustering_distance_cols = sampleDists,
         clustering_distance_rows = sampleDists)

# Here, we perform replicate correlation analysis using Pearson Correlation:

arld <- assay(rld)
logcpm <- cpm(arld[,order(colnames(arld))]) #Creates a matrix of normalized log2 counts per million values from the vsd object
logcpm_df <- as.data.frame(logcpm)
pcor <- rcorr(logcpm)
melted_pcor <- melt(pcor$r)
ggplot(data = melted_pcor, aes(x = Var1, y = Var2, fill = value)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ggtitle('Pearson Correlation')+geom_tile() + 
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.995, limit = c(0.99,1), space = "Lab", 
                       name = "Pearson\nCorrelation\nValue")

# Here, let's specify some genes of interest for highlighting in further analysis/plotting:

#genes_of_interest <- c('whatever you want')
gprin_genes <- c('Crym', 'Gprin3', 'Slco2a1')
colgalt_genes <- c('Colgalt2', 'Fezf2', 'Pcp4', 'Npnt')
L5b_genes <- c('Nefh', 'Bcl11b', 'Nrip3')
striat_genes <- c('Drd1', 'Drd2')
mito_genes <- c('Cox6c', 'Ndufa4', 'Cox6a1', 'Ndufb3', 'Uqcrh', 'Atp5g3')
all_genes <- c(gprin_genes, colgalt_genes, L5b_genes, striat_genes, mito_genes)

# Here, we make a heatmap that clusters our samples by expression of genes of interest from above.
# If you want to specify a new set of these genes, do that below:

pheatmap(subset(norm_counts[,order(colnames(norm_counts))], rownames(norm_counts) %in% c(genes_of_interest)),
         scale = 'row',
         color = colorRampPalette(c('black', 'mediumvioletred', 'darkorange', 'cornsilk'), space='rgb')(44), 
         cluster_cols = FALSE,
         main = 'Marker genes')

# Now, we finally run DE analysis
# We tell DESeq to compare an experimental sample against a control from our "sample_type"s
# Note that the control sample is listed last in the comparison:

# For determining the DE genes between Colgalt2 and Gprin3:

dds$group <- factor(paste0(dds$cell, dds$when, dds$type, dds$disease, dds$age, dds$ip_input))
design(dds) <- ~ group
dds <- DESeq(dds)
#head(resultsNames(dds))

getResults <- function(group, sample1, sample2) {
  res <- results(dds, contrast=c(group, sample1, sample2))
  res_drop <- as.data.frame(res)
  res_drop[is.na(res_drop)] <- 1
  return(res_drop)
}

res_whole_M1_gprin <- getResults('group', 'Gprin3M1NormWTPreIP', 'Gprin3WholeNormWTPreIP')

res_whole_M1_colgalt <- getResults('group', 'Colgalt2M1NormWTPreIP', 'Colgalt2WholeNormWTPreIP')

res_M1_gprin_colgalt <- getResults('group', 'Gprin3M1NormWTPreIP', 'Colgalt2M1NormWTPreIP')

res_whole_gprin_colgalt <- getResults('group', 'Gprin3WholeNormWTPreIP', 'Colgalt2WholeNormWTPreIP')

res_M1_gprin_IP_input <- getResults('group', 'Gprin3M1NormWTPreIP', 'Gprin3M1NormWTPreInput')

res_M1_colgalt_IP_input <- getResults('group', 'Colgalt2M1NormWTPreIP', 'Colgalt2M1NormWTPreInput')

res_M1_gprin_whole_colgalt <- getResults('group', 'Gprin3M1NormWTPreIP', 'Colgalt2WholeNormWTPreIP')

res_gprin_pre_als <- getResults('group', 'Gprin3WholeALSSODPreIP', 'Gprin3WholeALSWTPreIP')

res_gprin_post_als <- getResults('group', 'Gprin3WholeALSSODPostIP', 'Gprin3WholeALSWTPostIP')

res_colgalt_pre_als <- getResults('group', 'Colgalt2WholeALSSODPreIP', 'Colgalt2WholeALSWTPreIP')

res_colgalt_post_als <- getResults('group', 'Colgalt2WholeALSSODPostIP', 'Colgalt2WholeALSWTPostIP')

# MA plots:

plot_my_MA <- function(name_of_res, p_to_use = 'padj', gene_list = NULL) {
  of_interest <- name_of_res
  of_interest$isSig <- of_interest[,c(p_to_use)] < 0.05
  DESeq2::plotMA(of_interest[,c('baseMean', 'log2FoldChange', 'isSig')], cex = 0.6, 
                 main = deparse(substitute(name_of_res)))
  if (length(gene_list) > 1) {
    for (i in 1:length(gene_list)) {
      points(of_interest[gene_list[[i]],], col = i+2, pch = 16, cex = 0.6)
    } 
  } else {
    points(of_interest[gene_list,], col = 'yellow', pch = 16, cex = 0.6)
  }
  
  #text(of_interest[gene_list,], labels = mito_genes, col = 'green', pos = 4, cex = 0.6, font = 2)
}

plot_my_bar <- function(name_of_res, p_to_use = 'padj', gene_list = NULL) {
  of_interest <- name_of_res
  of_interest$isSig <- of_interest[,c(p_to_use)] < 0.05
  barplot(of_interest[c(gene_list),]$log2FoldChange,
          col = of_interest$isSig,
          names.arg = gene_list,
          cex.names = 0.6,
          las = 2)
}

# Now we plot the MA plot from these results,
# Coloring the genes with an adjusted pvalue less than 0.05 pink:

plot_my_MA(res_M1_gprin_colgalt, p_to_use='padj', gene_list = genes_of_interest)

plot_my_bar(res_gprin_post_als, 'pvalue', apop_genes)

# For external GO analysis, copy gene list to clipboard:

oi <- res_M1_gprin_colgalt
write_clip(rownames(oi[(oi$pvalue < 0.1 & oi$log2FoldChange < 0),]))

# Now we write the txt file with our results. We'll add the individual log counts for each replicate:

ind_counts <- as.data.frame(assay(rld))
relevant_cols <- ind_counts[,grep('IP', colnames(ind_counts))]
all_data <- merge(of_interest, relevant_cols, by ='row.names')

write.table(all_data, file = 'output.txt', sep ='\t', row.names = FALSE, col.names = TRUE)

# The End!

ddss <- DESeqDataSetFromMatrix(counts(dds,normalized=FALSE),
                               colData=newCol,
                               design = ~0+MouseLine+IP+MouseLine*IP)
ddss <- DESeq(ddss)

resultsNames(ddss)
res <- results(ddss,contrast = list(c("MouseLineColgalt2","MouseLineGprin3"))
               
               
