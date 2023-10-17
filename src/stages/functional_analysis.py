"uniprot.py"

import os
import pandas as pd

from src.common import http_post, http_request

##
# Depends upon transync and signet
##

WEBSITE_API = "https://rest.uniprot.org/"

def get_jobs_ids(genes_ids):
    # for id in genes_ids:
    return [
        http_post(f"{WEBSITE_API}/idmapping/run", data={"from": "Gene_Name","to": "UniProtKB","ids": gene_id}).json()['jobId']
        for gene_id in genes_ids
        ]  

def fetch_uniprots(job_ids):
    uni=[]
    for j in job_ids:
        r = http_request(f"{WEBSITE_API}/idmapping/status/{j}")
        data = r.json()
        for results in data['results']:
            try:
                if results['to']['organism']['scientificName'] == "Homo sapiens":
                    uni.append(([results['from'], results['to']['entryType'], results['to']['primaryAccession'], results['to']['comments'][0]['texts'][0]['value']]))
            except:
                print("An exception occurred")

    return pd.DataFrame(uni, columns=['Gene','Entry_type','Primary_acc','Function'])

def save_uniprot_data(target, uniprots_dataframe):
    os.mkdir(target+"/Uniprot")
    uniprots_dataframe.to_csv(target+"/Uniprot/Uniprot_transsynw_analysis.csv", index= False)


def start_uniprot(artefacts_path):
    print("RUNNING start_uniprot with params", artefacts_path)
    transync_genes_file = artefacts_path + "/Trrust_Analysis/transync_genes.csv"
    transsynw_genes = pd.read_csv(transync_genes_file)['Gene']
    jobs_ids = (get_jobs_ids(transsynw_genes))
    save_uniprot_data(artefacts_path,fetch_uniprots(jobs_ids))

# def get_genes_file(artefacts_path):
#     genes_file = artefacts_path + "/Trrust_Analysis/transync_genes.csv"
#     genes = pd.read_csv(genes_file)['Gene']
#     return genes
    
# def uniprot_search(artefacts_path,genes):
#     # genes_file = artefacts_path + "/Trrust_Analysis/transync_genes.csv"
#     # genes = pd.read_csv(genes_file)['Gene']
#     jobs_ids = (get_jobs_ids(genes))
#     save_uniprot_data(artefacts_path,fetch_uniprots(jobs_ids))


