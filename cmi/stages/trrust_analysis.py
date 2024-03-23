"trrust_analysis.py"
import pandas as pd
import numpy as np
import os
from cmi.common import clean_up
import logging


###
# Depends upon TransSyncW and SIGNET
###

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
file_handler = logging.FileHandler('logs.log')
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)


def read_trrust_db(db_file):
    trrust_data = pd.read_table(db_file ,header=None)
    trrust_data.columns =['Gene', 'Target', 'Action', 'Reference']
    return trrust_data.drop(['Reference'],axis = 1)

def combine_transync_outputs(cores_file, marker_file, target_genes_file):
    transync = pd.read_table(cores_file)
    transync.rename(columns = {'core':'Source'}, inplace = True)
    transynm = pd.read_table(marker_file)
    transynm.insert(3, 'Source', 'marker')
    transync_genes = pd.concat([transync, transynm])
    transync_genes.to_csv(target_genes_file)
    return transync_genes

def reform_signet_output(signet_out_file, signet_unique_genes_file):
    signet = pd.read_csv(signet_out_file)
    tf = signet['V1']
    tg = signet['V2']
    signet_data = pd.concat([tf,tg])
    signet_unique_genes = pd.DataFrame((signet_data.unique()),columns=['Gene'])
    signet_unique_genes['Source']="SIGNET"
    signet_unique_genes.to_csv(signet_unique_genes_file,index=False)
    return pd.DataFrame(signet_data,columns=['Gene'])

def store_nohits(target_genes_file, signet_unique_genes_file):
    transync_genes = pd.read_csv(target_genes_file)
    transync_genes = transync_genes[['Gene','Source']]
    signet_unique_genes = pd.read_csv(signet_unique_genes_file)
    combine=[transync_genes,signet_unique_genes]
    combined_data = pd.concat(combine)
    unique_combined_data = combined_data.drop_duplicates()
    
    return unique_combined_data


def analyse(trust_db, transync_combined, signet_reformed,unique_combined_data, artefacts_path): 
    trrust_transsynw_gene_match=trust_db[trust_db['Gene'].isin(transync_combined['Gene'])]
    trrust_transsynw_gene_match['Source']='TranSyn'
    trrust_transsynw_target_match=trust_db[trust_db['Target'].isin(transync_combined['Gene'])]
    trrust_transsynw_target_match['Source']='TranSyn'

    trrust_signet_gene_match=trust_db[trust_db['Gene'].isin(signet_reformed['Gene'])]
    trrust_signet_gene_match['Source']='SIGNET'
    trrust_signet_target_match=trust_db[trust_db['Target'].isin(signet_reformed['Gene'])]
    trrust_signet_target_match['Source']='SIGNET'
    
    trrust_analysis = [trrust_transsynw_gene_match,trrust_transsynw_target_match,trrust_signet_gene_match,trrust_signet_target_match]
    tdata = pd.concat(trrust_analysis)
    omitted = unique_combined_data[ ~unique_combined_data['Gene'].isin(tdata['Gene']) ]
    all_genes = [tdata,omitted]
    all_analysed_genes = pd.concat(all_genes)
    all_analysed_genes.fillna('No_data', inplace=True)
    all_analysed_genes.to_csv(artefacts_path + "/Trrust_Analysis/trrust_analysis.csv", index=False)
    return all_analysed_genes

def trrust_analysis(trust_db_file, artefacts_path):
    clean_up(artefacts_path+"/Trrust_Analysis")
    os.mkdir(artefacts_path+"/Trrust_Analysis")
    logger.info(f"RUNNING TRRUST_analysis with using: {trust_db_file, artefacts_path}")
    
    cores_file= artefacts_path+"/TransSynW/cores.tsv"
    markers_file = artefacts_path+"/TransSynW/markers.tsv"
    signet_file = artefacts_path+"/Signet/copaired2.csv"
    target_genes_file = artefacts_path+"/Trrust_Analysis/transsynw_genes.csv"
    signet_unique_genes_file = artefacts_path+"/Trrust_Analysis/signet_genes.csv"
    

    return analyse(read_trrust_db(trust_db_file),
    combine_transync_outputs(cores_file, markers_file, target_genes_file),
    reform_signet_output(signet_file,signet_unique_genes_file),
    store_nohits(target_genes_file, signet_unique_genes_file), artefacts_path)