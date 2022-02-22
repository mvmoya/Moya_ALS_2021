##########################################################################################################
#                                                                                                        #
#   This is the script used to load and analyze the files from the May-June 2019 Gprin3 and Colgalt2     #
#   IP experiments. It has been optimized to work for DESeq2 input and output.                           #
#                                                                                                        #
##########################################################################################################

# Here, we load the libraries that we'll need to run all of the analysis. Some of these can be deprecated

library(shiny)
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

##########################################################################################################
#                                              Load Functions                                            #
##########################################################################################################


getResults <- function(dds_object, group, sample2, sample1) {
  
  res <- results(dds_object, contrast=c(group, sample2, sample1))
  res_drop <- as.data.frame(res)
  res_drop[is.na(res_drop)] <- 1
  return(res_drop)
  
}

plot_my_MA <- function(name_of_res, p_to_use = 'padj', gene_list = NULL) {
  
  of_interest <- name_of_res
  of_interest$isSig <- of_interest[,c(p_to_use)] < 0.05
  
  DESeq2::plotMA(of_interest[,c('baseMean', 'log2FoldChange', 'isSig')], cex = 0.6, 
                 main = deparse(substitute(name_of_res)))
  
  if (length(gene_list) > 1) {
    
    colors <- vector()
    
    for (i in 1:length(gene_list)) {
      
      points(of_interest[gene_list[[i]],], col = i+2, pch = 16, cex = 0.6)
      colors <- c(colors, i+2)
    
    }
    
  } 
  
  else {
    
    points(of_interest[gene_list,], col = 'yellow', pch = 16, cex = 0.6)
  
  }
  
  legend_text <- sub("(list\\()", "", deparse(substitute(gene_list))) %>%
    sub("\\)", "", .) %>%
    strsplit(., ",")
  
  legend("bottomright", legend=unlist(legend_text), fill=colors, cex=0.6)

}

plot_my_bar <- function(name_of_res, p_to_use = 'padj', gene_list = NULL) {
  
  of_interest <- name_of_res
  of_interest$isSig <- of_interest[,p_to_use] < 0.05
  
  barplot(of_interest[c(gene_list),]$log2FoldChange,
          col = of_interest[c(gene_list),]$isSig,
          names.arg = gene_list,
          cex.names = 0.6,
          las = 2)
}

##########################################################################################################
#                                        Load/Structure Counts Files                                     #
##########################################################################################################

# Here, we set the working directory to wherever the counts files are all housed:


setwd('YOUR WORKING DIRECTORY')
current_wd <- getwd()

# This line gets a list of the files in the directory:

sampleFiles <- list.files(current_wd))

# Here, we create a list of sample conditions:

sampleCondition <- sub('.txt', '', sampleFiles)

# Now we make a reference table for the sample names from above and the corresponding filename:

sampleTable <- data.frame(fileName = sampleDir, condition = sampleCondition)

# Here we read the first file in the list and seed a data frame from it. I've changed the header
# names to be easier to reference, and have removed rows where gene names are NA values:

data.num <- read.delim(toString(sampleTable[1, c('fileName')]), header = TRUE)
data.num <- data.num[,2:3]
data.num <- data.num[complete.cases(data.num),]
colnames(data.num) <- c('genes', sampleCondition[1])

# Now, we add the remaining sample files to this table by merging the tables by gene name.
# This part loops through every filename specified in our reference table above:

for (i in 2:nrow(sampleTable)) {
  
  new_cols <- read.delim(toString(sampleTable[i,"fileName"]), header = TRUE)
  new_cols <- new_cols[,2:3]
  colnames(new_cols) <- c('genes', sampleCondition[i])
  data.num <- merge(data.num, new_cols, by = "genes")
  
}

# Now we make the rownames for the data frame into the gene names:

rownames(data.num) <- data.num[,1]
data.num[,1] <- NULL

##########################################################################################################
#                                         Create the DE Object                                           #
##########################################################################################################

# Here, we create the experimental information table. The user has to specify in "samplesDesired' which
# samples they want to include in making the DE Object:

sampleGroups <- unique(str_match(sampleCondition, '(.*)_\\d')[,2])
sampleGroups # Print out of the sampleGroups so you can copy and paste the names below:

samplesDesired <- c('ALS_Gprin3_SOD_TRAP_IP', 'ALS_Gprin3_WT_TRAP_IP')
columnsDesired <- data.num[,grep(paste(samplesDesired, collapse='|'), colnames(data.num))]

experiment <- sub('(Baseline|ALS).*', '\\1', colnames(columnsDesired))
cell <- sub('.*(Gprin3|Colgalt2).*', '\\1', colnames(columnsDesired))
disease <- sub(".*(WT|SOD).*", "\\1", colnames(columnsDesired))
ip_input <- sub(".*(IP|Input).*", "\\1", colnames(columnsDesired))

coldata <- data.frame(disease) #, cell, ip_input, age, region)
colnames(coldata) <- factor(c('disease')) #, 'cell', ip_input', 'age', 'region'))

# Now, we create the DESeq object, specifying that we want sample_type to be what groups our samples:

ddsPre <- DESeqDataSetFromMatrix(countData = columnsDesired, colData = coldata, design = ~ disease)
dds <- DESeq(ddsPre)

##########################################################################################################
#                                           Make Descriptive Plots                                       #
##########################################################################################################

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
color_scheme <- factor(as.numeric(factor(disease)))#*(as.numeric(factor(disease))+1)

rv <- rowVars(assay(vsd_of_interest))
select <- order(rv, decreasing=TRUE)#[seq_len(500)]
pca <- prcomp(t(assay(vsd_of_interest)[select,]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = disease, name = colnames(columnsDesired))
ggplot(data = d, aes_string(x = 'PC1', y = 'PC2', color = color_scheme, label = 'name')) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = name), size = 2, box.padding = 0.15, 
                  point.padding = 0.5, segment.color = 'grey50') +
  #geom_label_repel(aes(label = name), size = 2, box.padding = 0.35, 
  #                 point.padding = 0.5, segment.color = 'grey50') +
  xlab(paste('PC1: ', round(percentVar[1] * 100), '% variance')) +
  ylab(paste('PC2: ', round(percentVar[2] * 100), '% variance')) +
  coord_fixed() +
  ggtitle('PCA of samples')

View(pca$rotation[,1:2])

# Here, we perform heirarchical clustering of our samples based on gene expression:

sampleDists <- dist(t(assay(vsd_of_interest)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd_of_interest)
colnames(sampleDistMatrix) <- colnames(vsd_of_interest)
pheatmap(sampleDistMatrix, clustering_distance_cols = sampleDists,
         clustering_distance_rows = sampleDists)

# Here, we perform replicate correlation analysis using Pearson Correlation:

arld <- assay(rld)
logcpm <- cpm(arld[,order(colnames(sampleDistMatrix))]) #Creates a matrix of normalized log2 counts per million values from the vsd object
logcpm_df <- as.data.frame(logcpm)
pcor <- rcorr(logcpm)
melted_pcor <- melt(pcor$r)
ggplot(data = melted_pcor, aes(x = Var1, y = Var2, fill = value)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ggtitle('Pearson Correlation')+geom_tile() + 
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.99, limit = c(0.98,1), space = "Lab", 
                       name = "Pearson\nCorrelation\nValue")

# Here, let's specify some genes of interest for highlighting in further analysis/plotting:

# genes_of_interest <- c('whatever you want')
gprin_genes <- c('Crym', 'Gprin3', 'Slco2a1')
colgalt_genes <- c('Colgalt2', 'Fezf2', 'Pcp4', 'Npnt')
L5b_genes <- c('Nefh', 'Bcl11b', 'Nrip3')
striat_genes <- c('Drd1', 'Drd2')
mito_genes <- c('Cox6c', 'Ndufa4', 'Cox6a1', 'Ndufb3', 'Uqcrh', 'Atp5g3')
all_genes <- c(gprin_genes, colgalt_genes, L5b_genes, striat_genes, mito_genes)

candidate_gprin <- c('Lypd1', 'Ust')
candidate_colgalt <- c('Tmem150c', 'Vat1l')
candidate_both <- c('Parm1')


# Here, we make a heatmap that clusters our samples by expression of genes of interest from above.
# If you want to specify a new set of these genes, do that below after '%in%':

genes_to_plot <- c(candidate_gprin, candidate_colgalt, glial_genes, mito_genes)

pheatmap(subset(norm_counts[,order(colnames(norm_counts))], rownames(norm_counts) %in% genes_to_plot),
         scale = 'row',
         color = colorRampPalette(c('black', 'mediumvioletred', 'darkorange', 'cornsilk'), space='rgb')(44), 
         cluster_cols = FALSE,
         main = 'Marker genes')

##########################################################################################################
#                                           Perform DE Analysis                                          #
##########################################################################################################

# Now, we finally run DE analysis
# We tell DESeq to compare an experimental sample against a control from our "sample_type"s
# Note that the control sample is listed last in the comparison:

# For determining the DE genes between SOD and WT Gprin3:

res_M1_gprin_ALS <- getResults(dds, 'disease', 'SOD', 'WT') # control is second

##########################################################################################################
#                                              Plot DE Data                                              #
##########################################################################################################

# Now we plot the MA plot from these results, using the function above:

plot_my_MA(res_M1_gprin_ALS, p_to_use='padj', gene_list = list(mito_genes))

# And a bar graph of fold changes for specific genes:

plot_my_bar(res_M1_gprin_ALS, p_to_use='pvalue', mito_genes)

##########################################################################################################
#                                   Copy Gene Lists to Clipboard for GO                                  #
##########################################################################################################

# For external GO analysis, copy gene list to clipboard (be sure to adjust the comparisons as needed:

oi <- rownames(res_M1_gprin_ALS[(res_M1_gprin_ALS$pvalue < 0.05 & res_M1_gprin_ALS$baseMean > 100 & res_M1_gprin_ALS$log2FoldChange < -0.5),])
write_clip(oi)

##########################################################################################################
#                                           Write Output Table                                           #
##########################################################################################################

# Now we write the txt file with our results. We'll add the individual log counts for each replicate:

ind_counts <- as.data.frame(assay(rld))
relevant_cols <- ind_counts[,grep('IP', colnames(ind_counts))]
all_data <- merge(res_M1_gprin_ALS, relevant_cols, by ='row.names')

# Be sure to change the name of the file accordingly:

write.table(all_data, file = 'Gpring_SODvsWT_DESeq.txt', sep ='\t', row.names = FALSE, col.names = TRUE)

##########################################################################################################
#                                                 The End                                                #
##########################################################################################################
