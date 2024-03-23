"uniprot.py"

import os
import logging
import pandas as pd
from cmi.common import clean_up
from cmi.common import http_post, http_request


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
file_handler = logging.FileHandler('logs.log')
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)

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
                logger.info("An exception occurred")

    return pd.DataFrame(uni, columns=['Gene','Entry_type','Primary_acc','Function'])



def uniprot(artefacts_path, file):
    logger.info(f"Functional analysis started using: {artefacts_path, file}")
    file_path = artefacts_path + file
    file_name = str(file)
    file_name = file_name.split("/")
    file_name =file_name[2] or file

    genes = pd.read_csv(file_path)['Gene']
    jobs_ids = (get_jobs_ids(genes))
    hit = fetch_uniprots(jobs_ids)
    hit.to_csv(artefacts_path+"/Uniprot/uniprot_"+file_name,index= False)

def functional_analysis(artefacts_path, transsynw_genes_file, signet_genes_file):
    clean_up(artefacts_path+"/Uniprot/")
    os.makedirs(artefacts_path+"/Uniprot/")
    uniprot(artefacts_path, transsynw_genes_file)
    uniprot(artefacts_path, signet_genes_file)
    
    



