<<<<<<< HEAD
# Pipeline
=======
>>>>>>> 3cfc211 (created craft framework)
artefacts: Directory provided by user, where all the results would be saved.
stage: Part of the tool user wants to run.
params: Input arguments required by the stage.


*Commands*
- stage: all (TransSynW + PAGA + SIGNET + TRRUST + Cytoscape + Uniprot)
<<<<<<< HEAD

Eg: python main.py run all --artefacts ./artefacts/temp2 --params start.txt start_terminal.txt annotation.txt HPROGFPM HNES ./Signet/terminal.csv ./trrust_rawdata_human.tsv

- stage: mechanistic insights (TransSynW + PAGA)

Eg: python main.py run mechanistic_insights --artefacts ./artefacts/temp2 --params start.txt start_terminal.txt annotation.txt HPROGFPM HNES

- stage: grn inference (SIGNET)

Eg: python main.py run grn_inference --artefacts ./artefacts/temp --params ./Signet/terminal.csv

- stage: functional analysis (Uniprot)

Eg: python main.py run functional_analysis --artefacts ./artefacts/Results-2023-10-03_12-39-19

- stage: gene network (Cytoscape)

Eg: python main.py run generate_network --artefacts ./artefacts/Results-2023-10-03_12-39-19 --params ./artefacts/Results-2023-10-03_12-39-19/trrust_analysis.csv

- stage: trrust analysis (TRRUST)

Eg: python main.py run trrust_analysis --artefacts ./artefacts/Results-2023-10-03_12-39-19 --params ./trrust_rawdata_human.tsv 

- stage: generate_hypothesis (TransSynW)

=======
Eg: python main.py run all --artefacts ./artefacts/temp2 --params start.txt start_terminal.txt annotation.txt HPROGFPM HNES ./Signet/terminal.csv ./trrust_rawdata_human.tsv

- stage: mechanistic insights (TransSynW + PAGA)
Eg: python main.py run mechanistic_insights --artefacts ./artefacts/temp2 --params start.txt start_terminal.txt annotation.txt HPROGFPM HNES

- stage: grn inference (SIGNET)
Eg: python main.py run grn_inference --artefacts ./artefacts/temp --params ./Signet/terminal.csv

- stage: functional analysis (Uniprot)
Eg: python main.py run functional_analysis --artefacts ./artefacts/Results-2023-10-03_12-39-19
python main.py run functional_analysis --artefacts ./artefacts/stage_all_hffxprog transync_genes.csv signet_unique_gene_list.csv

- stage: gene network (Cytoscape)
Eg: python main.py run generate_network --artefacts ./artefacts/Results-2023-10-03_12-39-19 --params ./artefacts/Results-2023-10-03_12-39-19/trrust_analysis.csv
Eg: python main.py run generate_network --artefacts ./artefacts/Results-2023-10-03_12-39-19 --params ./artefacts/Results-2023-10-03_12-39-19/Trrust_Analysis/trrust_analysis.csv 

-stage: trrust analysis (TRRUST)
Eg: python main.py run trrust_analysis --artefacts ./artefacts/Results-2023-10-03_12-39-19 --params ./trrust_rawdata_human.tsv 

-stage: generate_hypothesis (TransSynW)
>>>>>>> 3cfc211 (created craft framework)
Eg: python main.py run generate_hypothesis --artefacts ./artefacts/temp2 --params start.txt start_terminal.txt annotation.txt HPROGFPM HNES
python main.py run generate_hypothesis --artefacts ${PWD}/artefacts/temp2 --params start.txt start_terminal.txt annotation.txt HPROGFPM HNES

setup file: downloads transsynw repo with code directory with the necessary codes
cannot chnage the name as for transsynw to run the path cannot be changed. Adding signet.py also in the code repository under Signet directory so that the dependencies are in one single folder.
