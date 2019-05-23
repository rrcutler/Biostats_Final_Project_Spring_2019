setwd('/Volumes/MyBook_3/BD_aging_project/RNAseq/All_tissues_analysis/DEseq2_runs/')
library(DESeq2)

# input count matrices are in "Figure1_General/RNA_Mapping_counting_DE/count_matrices"

####################################    Liver    #################################### 
# read in subread count matrix
my.liver1 <- read.table('count_matrices/Aging_Liver_counts_genes.txt',skip=1,header=T,sep="\t",stringsAsFactors=F)
my.liver <- my.liver1[,c(1,6:15)]
rownames(my.liver) <- my.liver[,1]

####################################   Heart   #################################### 
# read in subread count matrix
my.heart1 <- read.table('count_matrices/Aging_Heart_counts_genes.txt',skip=1,header=T,sep="\t",stringsAsFactors=F)
my.heart <- my.heart1[,c(1,6:15)]
rownames(my.heart) <- my.heart[,1]

#################################### Cerebellum #################################### 
# read in subread count matrix
# there were 2 nextseq runs based on poor clustering on flow cell
# will sum up count matrices
my.cereb1 <- read.table('count_matrices/Aging_cerebellum_counts_genes.txt',skip=1,header=T,sep="\t")
my.cereb2 <- read.table('count_matrices/Aging_cerebellum_v2_counts_genes.txt',skip=1,header=T,sep="\t")

my.cereb <- my.cereb1[,c(1,6:15)]
my.cereb[,3:11] <- my.cereb[,3:11] + my.cereb2[,7:15]
rownames(my.cereb) <- my.cereb[,1]

#################################### Olfactory Bulb #################################
# one of the 12mths samples was not analyzed
# read in subread count matrix
my.ob1 <- read.table('count_matrices/Aging_OlfactoryBulb_counts_genes.txt',skip=1,header=T,sep="\t",stringsAsFactors=F)
my.ob <- my.ob1[,c(1,6:14)]
rownames(my.ob) <- my.ob[,1]

####################################  NPCs pools  ###################################
# read in subread count matrix
my.npc1 <- read.table('count_matrices/Aging_NPCs_pool_counts_genes.txt',skip=1,header=T,sep="\t",stringsAsFactors=F)
my.npc <- my.npc1[,c(1,6:12)]
rownames(my.npc) <- my.npc[,1]


################
# create combined matrix
my.all <- cbind(my.liver,my.heart[,-c(1:2)],my.cereb[,-c(1:2)],my.ob[,-c(1:2)],my.npc[,-c(1:2)])
my.null <- which(apply(my.all[,3:length(my.all)], 1, sum) <= 1) # see deseq2 vignetter

# Now pull out the spike in genes
spikes.idx <- grep("ERCC-", rownames(my.all))
my.exclude <- union(my.null,spikes.idx)

my.filtered.matrix <- my.all[-my.exclude,3:length(my.all)]
rownames(my.filtered.matrix) <- my.all[-my.exclude,1]


# get age and cell type vectors
age <- as.numeric(c(rep(3,3),rep(12,3),rep(29,3) ,
                    rep(3,3),rep(12,3),rep(29,3),
                    rep(3,3),rep(12,3),rep(29,3),
                    rep(3,3),rep(12,2),rep(29,3),
                    rep(3,2),rep(12,2),rep(29,2))) # age in months

tissue <- c(rep("liver",9),rep("heart",9),rep("cereb",9),rep("OB",8),rep("NPCs",6))


# global design matrix
dataDesign = data.frame( row.names = colnames( my.filtered.matrix ), 
                         age = age, 
                         tissue = tissue )


#################################################################################################################
####             Normalize to get global age effect (age and tissues as independent covariates)               ###
#################################################################################################################
my.outprefix <- paste(Sys.Date(),"ALL_global_variance estimate","DESeq2_LINEAR_model_with_age",sep="_")


# get matrix using age as a modeling covariate
dds.1 <- DESeqDataSetFromMatrix(countData = my.filtered.matrix,
                              colData = dataDesign,
                              design = ~ age + tissue)

# run DESeq normalizations and export results
dds.deseq.1 <- DESeq(dds.1)

res.1 <- results(dds.deseq.1, name= "age") # added the name of the tested variable: doesn't seem to be taken correctly by default for FC

tissue.cts <- log2( counts(dds.deseq.1, normalize = TRUE) + 0.01)


# plot dispersion (variance modeled globally)
my.disp.out <- paste(my.outprefix,"_dispersion_plot.pdf")

pdf(my.disp.out)
plotDispEsts(dds.deseq.1)
dev.off()


# expression range
my.exp.out <- paste(my.outprefix,"_Normalized_counts_boxplot.pdf")

pdf(my.exp.out)
boxplot(tissue.cts,col=c(rep("coral",3),rep("blueviolet",3),rep("dodgerblue",3),
                         rep("coral",3),rep("blueviolet",3),rep("dodgerblue",3),
                         rep("coral",3),rep("blueviolet",3),rep("dodgerblue",3),
                         rep("coral",3),rep("blueviolet",2),rep("dodgerblue",3),
                         rep("coral",2),rep("blueviolet",2),rep("dodgerblue",2)),
        cex=0.5,ylab="Log2 DESeq2 Normalized counts", main = "ALL_GLOBAL",
        las=2, cex.axis = 0.5)  
dev.off()

### get the heatmap of aging changes at FDR5
## exclude NA
res.1 <- res.1[!is.na(res.1$padj),]

genes.aging <- rownames(res.1)[res.1$padj < 0.05]
my.num.aging <- length(genes.aging)

my.heatmap.out <- paste(my.outprefix,"_Heatmap_significant_genes.pdf")

pdf(my.heatmap.out)
my.heatmap.title <- paste("ALL"," aging singificant (FDR<5%), ",my.num.aging, " genes",sep="")
pheatmap(tissue.cts[genes.aging,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, cellwidth = 10)
dev.off()

# regress out the non age variance for plotting
full.model <- model.matrix(~ age + tissue, data = dataDesign) # all variables
fit <- lmFit(ExpressionSet(assayData=as.matrix(my.filtered.matrix)), full.model)
fit.eb <- eBayes(fit)
print(colnames(fit))
#[1] "(Intercept)" "age"         "tissueheart" "tissueliver" "tissueNPCs"  "tissueOB"   

### Regress out tissue ### 
#mod <- coefficients(fit)[,-c(1:2)] %*% t(fit$design[,-c(1:2)]) ### I keep only age (and intercept has cerebellum)
#
mod <- coefficients(fit)[,-2] %*% t(fit$design[,-2]) ### I keep only age (and intercept has cerebellum)
my.filtered.matrix.corrected <- my.filtered.matrix - mod

# do MDS analysis on tisue regressed data
mds.result <- cmdscale(1-cor(my.filtered.matrix.corrected,method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
x <- mds.result[, 1]
y <- mds.result[, 2]

my.colors <- c(rep("coral",3),rep("blueviolet",3),rep("dodgerblue",3),
               rep("coral",3),rep("blueviolet",3),rep("dodgerblue",3),
               rep("coral",3),rep("blueviolet",3),rep("dodgerblue",3),
               rep("coral",3),rep("blueviolet",2),rep("dodgerblue",3),
               rep("coral",2),rep("blueviolet",2),rep("dodgerblue",2))

my.mds.out <- paste(my.outprefix,"GLOBAL_aging_MDS_plot.pdf",sep="")

pdf(my.mds.out)
plot(x, y, xlab = "MDS dimension 1", ylab = "MDS dimension 2",main="Multi-dimensional Scaling",cex=2,pch=1)
points(x, y, pch=16,col=my.colors,cex=2)
legend("topleft",c("3m","12m","29m"),col=c("coral","blueviolet","dodgerblue"),pch=16,bty='n',pt.cex=1.5)
dev.off()


### get the heatmap of aging changes at FDR5 (tissue_regressed
## exclude NA

my.heatmap.out2 <- paste(my.outprefix,"_Heatmap_significant_genes_TISSUE_REGRESSED.pdf")

pdf(my.heatmap.out2, height = 10, width = 10, onefile=F)
my.heatmap.title <- paste("ALL TISSUE REGRESSED"," aging singificant (FDR<5%), ",my.num.aging, " genes",sep="")
pheatmap(my.filtered.matrix.corrected[genes.aging,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, cellwidth = 10)
dev.off()



# output result tables to files
my.out.ct.mat <- paste(my.outprefix,"_log2_counts_matrix.txt")
my.out.stats <- paste(my.outprefix,"_all_genes_statistics.txt")
my.out.fdr5 <- paste(my.outprefix,"_FDR5_genes_statistics.txt")
my.out.rdata <- paste(my.outprefix,"ALL_statistics.RData")

write.table(tissue.cts, file = my.out.ct.mat , sep = "\t" , row.names = T, quote=F)
write.table(res.1, file = my.out.stats , sep = "\t" , row.names = T, quote=F)
write.table(res.1[genes.aging,], file = my.out.fdr5, sep = "\t" , row.names = T, quote=F)

# save RData Object
my.aging_all.RNAseq.process <- res.1
save(res.1, file = my.out.rdata)