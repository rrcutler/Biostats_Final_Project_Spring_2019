setwd('/Volumes/BB_Backup_3//BD_aging_project/RNAseq/All_tissues_analysis/DEseq2_runs/Separate/')
source('RNAseq_analysis_functions_v4_forPCA.R')

# input count matrices are in "Figure1_General/RNA_Mapping_counting_DE/count_matrices"


####################################    Liver    #################################### 
# read in subread count matrix
my.liver1 <- read.table('/Volumes/BB_Backup_3/BD_aging_project/RNAseq/Liver/STAR/Aging_Liver_counts_genes.txt',skip=1,header=T,sep="\t",stringsAsFactors=F)
my.liver <- my.liver1[,c(1,6:15)]
rownames(my.liver) <- my.liver[,1]

# process RNAseq data and save RData object
my.liver.RNAseq.process <- process_aging_rnaseq("Liver", my.liver)
save(my.liver.RNAseq.process, file="RNA_seq_result_Liver_2018-03-09.RData")
##################################################################################### 



####################################   Heart   #################################### 
# read in subread count matrix
my.heart1 <- read.table('/Volumes/BB_Backup_3/BD_aging_project/RNAseq/Heart/STAR/Aging_Heart_counts_genes.txt',skip=1,header=T,sep="\t",stringsAsFactors=F)
my.heart <- my.heart1[,c(1,6:15)]
rownames(my.heart) <- my.heart[,1]

# process RNAseq data and save RData object
my.heart.RNAseq.process <- process_aging_rnaseq("Heart", my.heart)
save(my.heart.RNAseq.process, file="RNA_seq_result_Heart_2018-03-09.RData")
###################################################################################



#################################### Cerebellum #################################### 
# read in subread count matrix
# there were 2 nextseq runs based on poor clustering on flow cell
# will sum up count matrices
my.cereb1 <- read.table('/Volumes/BB_Backup_3/BD_aging_project/RNAseq/Cereb/1st_run/STAR/Aging_cerebellum_counts_genes.txt',skip=1,header=T,sep="\t")
my.cereb2 <- read.table('/Volumes/BB_Backup_3/BD_aging_project/RNAseq/Cereb/2nd_run/STAR/Aging_cerebellum_v2_counts_genes.txt',skip=1,header=T,sep="\t")

# process RNAseq data and save RData object
my.cereb.RNAseq.process <- process_aging_rnaseq("Cerebellum", my.cereb)
save(my.cereb.RNAseq.process, file="RNA_seq_result_cereb_2018-03-09.RData")
####################################################################################



#################################### Olfactory Bulb #################################
# one of the 12mths samples was not analyzed
# read in subread count matrix
my.ob1 <- read.table('/Volumes/BB_Backup_3/BD_aging_project/RNAseq/OB/STAR/Aging_OlfactoryBulb_counts_genes.txt',skip=1,header=T,sep="\t",stringsAsFactors=F)
my.ob <- my.ob1[,c(1,6:14)]
rownames(my.ob) <- my.ob[,1]

# process RNAseq data and save RData object
my.ob.RNAseq.process <- process_aging_rnaseq("OlfactoryBulb", my.ob, reps.3=3, reps.12=2, reps.29=3)
save(my.ob.RNAseq.process, file="RNA_seq_result_OB_2018-03-09.RData")
##################################################################################### 



####################################  NPCs pools  ###################################
# read in subread count matrix
my.npc1 <- read.table('/Volumes/BB_Backup_3/BD_aging_project/RNAseq/NPC_Pool/STAR/Aging_NPCs_pool_counts_genes.txt',skip=1,header=T,sep="\t",stringsAsFactors=F)
my.npc <- my.npc1[,c(1,6:12)]
rownames(my.npc) <- my.npc[,1]

# process RNAseq data and save RData npcject
my.npc.RNAseq.process <- process_aging_rnaseq("NPCs", my.npc, reps.3=2, reps.12=2, reps.29=2)
save(my.npc.RNAseq.process, file="RNA_seq_result_NPCs_2018-03-09.RData")
##################################################################################### 
