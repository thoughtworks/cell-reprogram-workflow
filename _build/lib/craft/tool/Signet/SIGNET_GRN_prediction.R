
getwd()
args= commandArgs(trailingOnly = TRUE) 
in_path = args[1]


# code_path =paste(wd+'/src/tool/Signet/')

library(reshape2)
library(umap)
library(pheatmap)
library(igraph)
library(GGally) # ggplot2 >= 3.3.0
library(ggplot2)
library(RcisTarget)
library(AUCell)
library(RcisTarget)

source('./craft/tool/Signet/SIGNET.R')



# Select motif database to use (annotations)
data("motifAnnotations_hgnc")

# Import the motif databases for RcisTarget
motifRankings <- importRankings('./craft/tool/Signet/hg19-tss-centered-10kb-10species.mc9nr.feather')
gene <- colnames(motifRankings)


ac_ntf <- read.table(paste(in_path,"/Signet/Intermediate/co_fc.txt",sep=""))
ac_tf <- read.table(paste(in_path,"/Signet/Intermediate/co_tf_fc.txt",sep=""))
gene_ntf <- read.csv(paste(in_path,"/Signet/Intermediate/gene_ntf.csv",sep=""))
gene_tf <- read.csv(paste(in_path,"/Signet/Intermediate/gene_tf.csv",sep=""))

# Merge data
coexpressed <- Merge(ac_ntf,ac_tf,gene_ntf,gene_tf,gene)
# gESD test and outliers screening
Outliers <- Screen(coexpressed)


# Transform to gene list
copaired <- Copaired(Outliers)
write.csv(copaired, (paste(in_path,"/Signet/copaired1.csv",sep="")))
genesets <- Genesets(copaired)
write.csv(genesets,(paste(in_path,"/Signet/genesets1.csv",sep="")))
genesets_list <- levels(as.factor(copaired$V1))
write.csv(genesets_list,(paste(in_path,"/Signet/genesets_list1.csv",sep="")))

# Motif enrichment analysis
# Three species selections are available: 'Homo sapiens','Mus musculus','Drosophila melanogaster'
# Choose right selection during analysis

Regulons <- Prune(Outliers, genesets, motifRankings, motifAnnotations, species = 'Homo sapiens')
copaired2 <- Copaired2(Regulons)
write.csv(copaired2,(paste(in_path,"/Signet/copaired2.csv",sep="")))

genesets2 <- Genesets(copaired2)
write.csv(genesets2,(paste(in_path,"/Signet/genesets2.csv",sep="")))

genesets_list2 <- levels(as.factor(copaired2$V1))
write.csv(genesets_list2,(paste(in_path,"/Signet/genesets_list2.csv",sep="")))


print("signet_R_done!")
