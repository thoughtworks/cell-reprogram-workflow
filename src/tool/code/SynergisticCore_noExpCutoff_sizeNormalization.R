#!/usr/bin/env Rscript

## Usage: Rscript code/SynergisticCore_noExpCutoff_sizeNormalization.R sample_data_Gokce2016.csv ::Oligo::Astro:: . sample_data_Gokce2016_cluster.csv hESC_E-MTAB-6819_H9_naive_smartseq2.Robj Mouse sample_data::sample_cluster

## Import all.txt files in the current directory
args = commandArgs(TRUE)
infile = args[1]  # data.csv
outDir = args[3]
clusterfile = args[4]  # cluster file
donorC = args[5]  # initial cell type
species = args[6]
orgname = args[7]
#setwd('src/tool/')
## Load necessary libraries
path = 'src/tool/'
source(paste(path,'code/SynergisticCore_library_noExpCutoff_pioneers_synergyBG.R',sep="")) # 5 + all 3-4-5 nonJSDs
source(paste(path,'code/scJSD.R',sep=""))
library(Rcpp)
sourceCpp(paste(path,'code/SynergisticCore.cpp',sep=""))
sourceCpp(paste(path,'code/JSD.cpp',sep=""))
require(gtools)
library(Matrix)
library(tibble)
library(dplyr)
library(stringr)
library(purrr)

target.classes = str_split(args[2], "::", simplify = TRUE) %>% as.vector
merge=target.classes[length(target.classes)]
target.classes = target.classes[-length(target.classes)][-1]
orgname = str_split(args[7], "::", simplify = TRUE) %>% as.vector

print(target.classes)
print(merge)
print (species)
print(orgname)

## Load start pop
if(! file.exists((paste0("startpop/",donorC)))){
  geneexp=orgname[1]
  clust=orgname[2]
  donor_org=orgname[3]
  inipop = read.delim(donorC,sep="\t",header=TRUE,stringsAsFactors=FALSE,row.names=1)
  inipop = log2(inipop+1)
  donorE = rowMeans(inipop)
  #save(donorE, file='initial_gene_expression.Robj')
} else {
  donor_org=donorC
  geneexp=orgname[1]
  clust=orgname[2]
  load((paste0("startpop/",donorC)))
}

## Load data
ftype = system(paste0("file -b --mime-type ",infile," | sed 's|/.*||'"),intern=TRUE)
if(ftype=='text'){
  M_all = read.delim(infile,sep="\t",header=TRUE,stringsAsFactors=FALSE,row.names = 1)
} else {
  M_all = readRDS(infile)
}
rownames(M_all) = toupper(rownames(M_all))


# TF subsetting

tfs = read.csv(paste0("src/tool/code/TF_",species,".txt"),header=FALSE,sep="\t",quote="")
M = M_all[toupper(rownames(M_all)) %in% toupper(tfs[,1]),]
print(dim(M))


## Assign cluster IDs
map = read.delim(clusterfile,sep="\t",header=FALSE,stringsAsFactors=FALSE)
print(map)

m = match(colnames(M), as.character(map[1,]))
print(m)
if(sum(is.na(m))==0){
  colnames(M) = map[2,m]
}else{
  stop("cluster-unassigned cells exist")
}


# merge subpopulations
if(merge=="yes"){
  colnames(M)[colnames(M) %in% target.classes] = paste0(target.classes,collapse="_")
  target.classes = paste0(target.classes,collapse="_")
}

## All classes
classes = unique(colnames(M));


## Identify most synergistic cores.
for(cla in target.classes){
  if(cla != 'undefined'){ 
    print(cla)
    computeMostSynergisticCores(data=M, classes=cla, infile=infile, donorC=donorC, outDir=outDir)
  }
}

## output
fil2 = paste0(outDir,"/analysis_summary.txt")
write.table(paste("Starting cell type:",donor_org, sep=" "), file=fil2,sep="\t",row.names=FALSE,quote=FALSE, col.names=FALSE,append=TRUE)
write.table(paste("Gene expression matrix:",geneexp, sep=" "), file=fil2,sep="\t",row.names=FALSE,quote=FALSE, col.names=FALSE,append=TRUE)
write.table(paste("Cluster file:",clust, sep=" "), file=fil2,sep="\t",row.names=FALSE,quote=FALSE, col.names=FALSE,append=TRUE)
write.table(paste("Species:",species, sep=" "), file=fil2,sep="\t",row.names=FALSE,quote=FALSE, col.names=FALSE,append=TRUE)
write.table(paste("Analysed subpopulations:",paste0(target.classes,collapse = ", "), sep=" "), file=fil2,sep="\t",row.names=FALSE,quote=FALSE, col.names=FALSE,append=TRUE)

core_files = list.files(path=outDir, pattern = "core",full.names = TRUE)

print(core_files)

cores = map_dfr(core_files, function(x) {
  read.delim(x, sep = ',', stringsAsFactors = F)
})
write.table(cores, paste0(outDir,"/cores.tsv"), sep="\t", row.names = F)
unlink(core_files)

## Markers
gens = read.csv('src/tool/code/all_marker_genes.txt',sep="\t",quote="",header=FALSE)
if(sum(toupper(rownames(M_all)) %in% toupper(gens[,1])) < 1000 ){stop("gene expression matrix does not contain enough marker genes")}

print(">Marker computation")
M_2 = M_all[toupper(rownames(M_all)) %in% toupper(gens[,1]), ]
print(dim(M_2))
M = M[rowSums(M_2)!=0,,drop=FALSE]



# sizeNormalization
#write.csv(M_2,"/Users/avani_mahadik/Downloads/scp_pancreatic_islets/early_data_transrun/before_size_normalisation.csv")
# print(M)
M_2 = sweep(M_2,2,colSums(M_2),`/`)
options("scipen"=9999, "digits"=10)
#write.csv(M_2,"/Users/avani_mahadik/Downloads/scp_pancreatic_islets/early_data_transrun/after_sweep.csv")



# print(M)

M_2[is.na(M_2)] = 0
M_2 = M_2[rowSums(M_2)!=0,,drop=FALSE]


JSD_rank_threshold = 10
tst = 0
for(cla in target.classes){
  if(cla != 'undefined'){
    tst = tst+1
    #print("Running C++ JSD...")
    print("Printing JSD multiple times.......")
    print(tst)
    JSD_value = scJSD(M_2, cla)
    #write.csv(JSD_value,"/Users/avani_mahadik/Documents/data_from_code/not_working/JSD_value.csv")
    print("Printing JSD matrix dim")
    print(dim(JSD_value))
    print("##### printing JSD_value ##### ")
    print(JSD_value)
    sJSD_expr_score = scJSD_score(JSD_value)
    print("##### printing sJSD_expr_score ######")
    print(sJSD_expr_score)
    #write.csv(sJSD_expr_score,"/Users/avani_mahadik/Documents/data_from_code/not_working/JSD_value_aftersJSD_expr.csv")
    temp = sJSD_expr_score[1:JSD_rank_threshold]
    temp = cbind('Gene'=names(temp),'scJSD'=round(temp,2))
    fil = paste0(outDir,"/markers_",cla,".txt")
    write.table(temp, file=fil,sep=",",row.names=FALSE,quote=FALSE, col.names=TRUE)
  }	
}

# output
marker_files = list.files(path=outDir, pattern = "markers", full.names = TRUE)
markers = map_dfr(marker_files, function(x) {
  temp_df = read.delim(x, sep = ',', stringsAsFactors = F)
  pop = str_match(x, "markers_(.+)\\.txt")[2]
  temp_df = temp_df %>% mutate(Subpopulation = pop)
})
write.table(markers, paste0(outDir,"/markers.tsv"), sep="\t", row.names = F)
unlink(marker_files)



