"generate_network.py"

import py4cytoscape as p4c
import pandas as pd

def unique_genes_list(trust_analysis_out_file):
    genes_list = pd.read_csv(trust_analysis_out_file)
    G = genes_list['Gene']
    T = genes_list['Target']
    merge = [G,T]
    merge = pd.concat(merge)
    merge.to_list()
    return merge.unique()
    

def generate(genes, trrust_analysis):
    p4c.cytoscape_ping()
    p4c.cytoscape_version_info()
    tnodes = pd.DataFrame(data={'id': genes})
    print(tnodes)
    tedges = pd.DataFrame(data={'source': trrust_analysis['Gene'], 'target': trrust_analysis['Target'], 'interaction': trrust_analysis['Action']})
    print(tedges)
    p4c.create_network_from_data_frames(tnodes, tedges, title="Network", collection="Trrust_Database")
    p4c.set_visual_style('Marquee')



def generate_network(trrust_out):
    print("RUNNING generate_network with params", trrust_out)

    generate(unique_genes_list(trrust_out), pd.read_csv(trrust_out))