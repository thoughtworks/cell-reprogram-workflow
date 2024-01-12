setup.sh file: User needs to run this file to start. This file downloads the dependencies and runs the pip installation of craft. 
artefacts: Directory provided by user, where all the results would be saved.
stage: Part of the tool user wants to run.
params: Input arguments required by the stage.

python version used: 3.9
Cytoscape version recommended: 3.9.1. Cytoscape needs to open in the background while running the workflow.
PathLinker app to be added to cytoscape using App manager in the software.
R requirements: 
libraries in rstudio: reshape2, umap, pheatmap, igraph ,GGally, ggplot2 >= 3.3.0 ,RcisTarget ,AUCell ,RcisTarget

All the input files should be saved in a folder called data.
Inputs files user needs to create:
Gene expression files
1. Starting cell population: create a .txt file with gene expression data of starting cell population. The rows should be gene names and columns should be sample IDs.
Eg: start.txt
2. Starting cell and terminal cell population combined: create a .txt file with gene expression data of both starting cell and terminal cell population. The rows should be gene names and columns should be sample IDs.
Eg: start_terminal.txt
![example start and terminal population combined](images/eg_start+ter_data_pic.png)
3. Terminal cell population: create a .csv file with gene expression data of starting cell population. The rows should be gene names and columns should be sample IDs.
Eg: terminal.csv
4. annotation.txt: create a .txt file with lable IDs of the samples from starting cell and terminal cell population combined expression data. First row being the sample IDS similar to the starting cell and terminal cell population combined file and second two should be the cluster IDs of the population. Also one of the sample ID and cluster ID should be matching.
![example annotation image](images/eg_annotation_pic.png)

The starting cell population and terminal cell population cluster IDs to be enerted as parameters should match the one in the annotations files.


**Commands**
<u>stage</u>: all (TransSynW + PAGA + SIGNET + TRRUST + Cytoscape + Uniprot)
#craft run all --artefacts ./artefacts/[directory_name] --params [start_cell population] [start and terminal_cell population] [annotation file] [terminal cell cluster ID] [startaing cell cluster ID] ./data/terminal.csv ./data/trrust_rawdata_human.tsv

Eg: craft run all --artefacts ./artefacts/temp --params start.txt start_terminal.txt annotation.txt HPROGFPM HNES ./data/terminal.csv ./data/trrust_rawdata_human.tsv

<u>stage</u>: generate_hypothesis (TransSynW)
# craft run all --artefacts ./artefacts/[directory_name] --params [start_cell population] [start and terminal_cell population] [annotation file] [terminal cell cluster ID]

Eg: craft run generate_hypothesis --artefacts ./artefacts/[directory_name] --params start.txt start_terminal.txt annotation.txt HPROGFPM

<u>stage</u>: mechanistic insights (TransSynW + PAGA)
# craft run all --artefacts ./artefacts/[directory_name] --params [start_cell population] [start and terminal_cell population] [annotation file] [terminal cell cluster ID] [startaing cell cluster ID]

Eg: craft run mechanistic_insights --artefacts ./artefacts/temp --params start.txt start_terminal.txt annotation.txt HPROGFPM HNES

<u>stage</u>: grn inference (SIGNET)
# craft run grn_inference --artefacts ./artefacts/[directory_name] --params ./data/terminal.csv
Eg: craft run grn_inference --artefacts ./artefacts/temp --params ./data/terminal.csv

<u>stage</u>: functional analysis (Uniprot)
# craft run functional_analysis --artefacts ./artefacts/[directory_name] transync_genes.csv signet_unique_gene_list.csv
Eg: craft run functional_analysis --artefacts ./artefacts/temp transync_genes.csv signet_unique_gene_list.csv

<u>stage</u>: gene network (Cytoscape)
# craft run generate_network --artefacts ./artefacts/[directory_name] --params ./artefacts/[directory_name]/Trrust_Analysis/trrust_analysis.csv

Eg: craft run generate_network --artefacts ./artefacts/temp --params ./artefacts/temp/Trrust_Analysis/trrust_analysis.csv 

<u>stage</u>: trrust analysis (TRRUST)
# craft run trrust_analysis --artefacts ./artefacts/[directory_name] --params ./data/trrust_rawdata_human.tsv
Eg: craft run trrust_analysis --artefacts ./artefacts/temp --params ./data/trrust_rawdata_human.tsv 
