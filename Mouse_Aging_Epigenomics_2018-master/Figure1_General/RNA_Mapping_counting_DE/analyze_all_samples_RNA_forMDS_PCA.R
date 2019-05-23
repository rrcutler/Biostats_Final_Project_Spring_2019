setwd('/Volumes/MyBook_3/BD_aging_project/RNAseq/All_tissues_analysis/')
library(DESeq2)

# global clustering: do norm cross all samples
# generate global MDS/PCA

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



my.all <- cbind(my.liver,my.heart[,-c(1:2)],my.cereb[,-c(1:2)],my.ob[,-c(1:2)],my.npc[,-c(1:2)])

my.null <- which(apply(my.all[,3:length(my.all)], 1, sum) <= 1) # see deseq2 vignette

# Now pull out the spike in genes
spikes.idx <- grep("ERCC-", rownames(my.all))
my.exclude <- union(my.null,spikes.idx)

my.filtered.matrix <- my.all[-my.exclude,3:length(my.all)]
rownames(my.filtered.matrix) <- my.all[-my.exclude,1]

age <- as.numeric(c(rep(3,3),rep(12,3),rep(29,3) ,
                    rep(3,3),rep(12,3),rep(29,3),
                    rep(3,3),rep(12,3),rep(29,3),
                    rep(3,3),rep(12,2),rep(29,3),
                    rep(3,2),rep(12,2),rep(29,2))) # age in months

tissue <- c(rep("liver",9),rep("heart",9),rep("cereb",9),rep("OB",8),rep("NPCs",6))
# design matrix
dataDesign = data.frame( row.names = colnames( my.filtered.matrix ), age = age, tissue = tissue )

# get matrix using age as a modeling covariate
dds <- DESeqDataSetFromMatrix(countData = my.filtered.matrix,
                              colData = dataDesign,
                              design = ~ age + tissue)


dds.deseq <- DESeq(dds)

# get normalized data
tissue.cts <- log2( counts(dds.deseq, normalize = TRUE) + 0.01)

mds.result <- cmdscale(1-cor(tissue.cts,method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
x <- mds.result[, 1]
y <- mds.result[, 2]

my.colors <- c(rep("coral",3), rep("blueviolet", 3),rep("dodgerblue",3),
               rep("coral",3), rep("blueviolet", 3),rep("dodgerblue",3),
               rep("coral",3), rep("blueviolet", 3),rep("dodgerblue",3),
               rep("coral",3), rep("blueviolet", 2),rep("dodgerblue",3),
               rep("coral",2), rep("blueviolet", 2),rep("dodgerblue",2))

# pch: by tissues
my.pchs <- c(rep(8,9),rep(14,9),rep(5,9),rep(1,8),rep(11,6))

# NPC - 11
# Heart - 14
# Cere - 5
# Liver - 8
# OB - 1

pdf("2017-04-07_MDS_RNAseq_DESeq_norm_together_BIGPTS.pdf")
plot(x, y, xlab = "MDS dimension 1", ylab = "MDS dimension 2",main="Multi-dimensional Scaling",cex=4,col=NULL)
points(x, y, pch=my.pchs,col=my.colors,cex=4)
legend("topleft",c("NPCs","Cerebellum","Olfactory bulb","Heart","Liver"),pch=c(11,5,1,14,8),col="grey",bty='n',pt.cex=1)
legend("topright",c("3m","12m","29m"),col=c("coral","blueviolet","dodgerblue"),pch=16,bty='n',pt.cex=1)
dev.off()


##### do PCA analysis
my.pos.var <- apply(tissue.cts,1,var) >0
my.pca <- prcomp(t(tissue.cts[my.pos.var,]),scale = TRUE)
x <- my.pca$x[,1]
y <- my.pca$x[,2]

my.summary <- summary(my.pca)

my.pca.out <- paste(Sys.Date(),"PCA_RNAseq_DESeq_norm_together_BIGPTS.pdf",sep="")

pdf(my.pca.out)
plot(x,y,pch = 16, cex=4, 
     xlab = paste('PC1 (', round(100*my.summary$importance[,1][2],1),"%)", sep=""),
     ylab = paste('PC2 (', round(100*my.summary$importance[,2][2],1),"%)", sep=""),
     cex.lab = 1, col = NULL) 
points(x,y, cex= 4, lwd = 1.5, col=my.colors, pch = my.pchs)
legend("bottomright",c("NPCs","Cerebellum","Olfactory bulb","Heart","Liver"),pch=c(11,5,1,14,8),col="grey",bty='n',pt.cex=1)
legend("bottomleft",c("3m","12m","29m"),col=c("coral","blueviolet","dodgerblue"),pch=16,bty='n',pt.cex=1)
dev.off()