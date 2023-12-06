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

def save_uniprot_data_transyn(target, uniprots_dataframe):
    os.makedirs(target+"/Uniprot/transync")
    uniprots_dataframe.to_csv(target+"/Uniprot/transync/Uniprot_transsynw_analysis.csv", index= False)

def save_uniprot_data_signet(target, uniprots_dataframe):
    os.makedirs(target+"/Uniprot/signet")
    uniprots_dataframe.to_csv(target+"/Uniprot/signet/Uniprot_signet_analysis.csv", index= False)



def start_uniprot_transyn(artefacts_path):
    print("RUNNING start_uniprot for transyn data with params", artefacts_path)
    transync_genes_file = artefacts_path + "/Trrust_Analysis/transync_genes.csv"
    transsynw_genes = pd.read_csv(transync_genes_file)['Gene']
    jobs_ids = (get_jobs_ids(transsynw_genes))
    save_uniprot_data_transyn(artefacts_path,fetch_uniprots(jobs_ids))

def start_uniprot_signet(artefacts_path):
    print("RUNNING start_uniprot for signet data with params", artefacts_path)
    signet_file = artefacts_path + "/Trrust_Analysis/signet_unique_gene_list.csv"
    signet_genes = pd.read_csv(signet_file)['Gene']
    jobs_ids = (get_jobs_ids(signet_genes))
    save_uniprot_data_signet(artefacts_path,fetch_uniprots(jobs_ids))



