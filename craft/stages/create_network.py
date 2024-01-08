"create_network.py"

import py4cytoscape as p4c
import pandas as pd
import os
import json
import requests
import logging
from craft.common import clean_up

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
file_handler = logging.FileHandler('logs.log')
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)

def unique_node_list(trust_analysis_out_file):
    trrust_analysis_data = pd.read_csv(trust_analysis_out_file)
    genes = trrust_analysis_data['Gene']
    targets = trrust_analysis_data['Target']
    node_data = [genes,targets]
    node_data = pd.concat(node_data)
    node_data.to_list()
    node_data = node_data.unique()
    return node_data

def generate_network(node_data, trrust_analysis_data, artefacts_path):
    clean_up(artefacts_path+"/Cytoscape/")
    os.mkdir(artefacts_path+"/Cytoscape/")
    nodes = pd.DataFrame(data={'id': node_data})
    edges = pd.DataFrame(data={'source': trrust_analysis_data['Gene'], 'target': trrust_analysis_data['Target'], 'interaction': trrust_analysis_data['Action']})
    p4c.create_network_from_data_frames(nodes, edges, title= "network_from_trrust_output" , collection="Analysis_using_cytoscape")
    p4c.export_network(artefacts_path+"/Cytoscape/"+'network_from_trrust_output', type='SIF',overwrite_file=True) 

def get_TFs_data(artefacts_path):
    signet = pd.read_csv(artefacts_path+"/Signet/copaired2.csv")
    cores = pd.read_table(artefacts_path+"/TransSynW/cores.tsv")
    signet_tf = signet['V1'].unique()
    signet_tf =' '.join(signet_tf)
    transsynw_core = cores['Gene'].unique()
    transsynw_core =' '.join(transsynw_core)
    TFs = transsynw_core + ' ' + signet_tf
    return str(TFs)
  
def get_targets_data(artefacts_path):
    signet = pd.read_csv(artefacts_path+"/Signet/copaired2.csv")
    markers = pd.read_table(artefacts_path+"/TransSynW/markers.tsv")
    signet_targets = signet['V2'].unique()
    signet_target = signet_targets.tolist()
    signet_target =' '.join(signet_target)
    transsynw_marker = markers['Gene'].unique()
    transsynw_marker =' '.join(transsynw_marker) 
    targets = transsynw_marker + ' ' + signet_target
    return str(targets)

def preprocess_sif(artefacts_path):
    path = artefacts_path +"/Cytoscape/"
    network_file = 'network_from_trrust_output.sif'
    p4c.import_network_from_file(file = path + network_file)
    p4c.networks.set_current_network(network_file) 
    p4c.select_nodes(['No_data'], by_col='name')
    p4c.delete_selected_nodes()
    preprocessed_file = p4c.export_network(path + 'preprocessed_' + network_file, type='SIF',overwrite_file=True)
    return(preprocessed_file)

def call_api(params,network_file):
    logger.debug(f'Sources (TFs): {params["sources"]}\nTargets (Marker + Target genes: {params["targets"]}')
    params["k"] = 100
    params["treatNetworkAsUndirected"] = False
    params["allowSourcesTargetsInPaths"] = True
    params["includeTiedPaths"] = True
    params["skipSubnetworkGeneration"] = False
    p4c.networks.set_current_network(network_file)
    p4c.clear_selection()
    headers = {'Content-Type': 'application/json', 'Accept': 'application/json'}
    url = "http://localhost:1234/pathlinker/v1/currentView/run"
    result_json = requests.request("POST", 
                          url,
                          data = json.dumps(params),
                          params = None,
                          headers = headers)


def get_regulators(TFs, targets, artefacts_path):
    path = artefacts_path +"/Cytoscape/"
    network_file = 'network_from_trrust_output.sif'
    p4c.networks.set_current_network(network_file)
    params = {}
    params['sources'] = TFs
    params['targets'] = targets
    call_api(params,network_file)
    p4c.export_network(filename =  artefacts_path+"/Cytoscape/" + 'PathLinker_TF_Targets_Paths.sif', type='SIF',overwrite_file=True)

def get_master_regulators(TFs,artefacts_path):
    path = artefacts_path +"/Cytoscape/"
    network_file = 'network_from_trrust_output.sif'
    p4c.networks.set_current_network(network_file) 
    TFs_list = TFs.split(' ')
    p4c.clear_selection()
    p4c.select_nodes(TFs_list, by_col='name')
    p4c.invert_node_selection()
    master_raw_list = p4c.get_selected_nodes()
    master_raw = ' '.join(master_raw_list)

    p4c.clear_selection()
    params = {}
    params['sources'] = master_raw
    params['targets'] = TFs
    call_api(params, network_file)
    p4c.export_network(filename = artefacts_path+"/Cytoscape/"+'PathLinker_MasterRegulators_TFs_Paths.sif', type='SIF',overwrite_file=True)

def merge_networks(artefacts_path):
    p4c.import_network_from_file(file = artefacts_path +"/Cytoscape/"+ 'PathLinker_TF_Targets_Paths.sif')
    p4c.set_visual_style(style_name='Directed', network='PathLinker_TF_Targets_Paths.sif')
    p4c.import_network_from_file(file = artefacts_path +"/Cytoscape/"+ 'PathLinker_MasterRegulators_TFs_Paths.sif')
    p4c.set_visual_style(style_name='Directed', network='PathLinker_MasterRegulators_TFs_Paths.sif')
    merged_network_file = 'MasterRegulators_TFs_Targets_Merged_Network.sif'

    try:
        p4c.tools.merge_networks(sources=['PathLinker_TF_Targets_Paths.sif', 'PathLinker_MasterRegulators_TFs_Paths.sif'], 
                                        title=merged_network_file, operation='union')
    except TypeError:
        logger.info(f'Error thrown by p4cytoscape during merging of networks. The merged network is created in Cytoscape though!')
    p4c.export_network(filename =  artefacts_path+ "/Cytoscape/MasterRegulators_TFs_Targets_Merged_Network" , type='SIF',overwrite_file=True)


def set_visual_style(merged_network_file, TFs, targets,  artefacts_path):
        
        p4c.import_network_from_file(file = artefacts_path +"/Cytoscape/"+ merged_network_file)
        p4c.networks.set_current_network(merged_network_file)
        p4c.set_visual_style(style_name='Directed', network=merged_network_file)
        all_genes_set = set(p4c.get_all_nodes())
        TFs_list = TFs.split(' ')
        TFs_in_network = list(set.intersection(all_genes_set, set(TFs_list)))
        p4c.style_bypasses.set_node_color_bypass(TFs_in_network, '#00FFFF')
        p4c.style_bypasses.set_node_shape_bypass(TFs_in_network, 'Diamond')
        targets_list = targets.split(' ')
        targets_in_network = list(set.intersection(all_genes_set,set(targets_list)))
        p4c.style_bypasses.set_node_color_bypass(targets_in_network, '#FFDF00')
        p4c.style_bypasses.set_node_shape_bypass(targets_in_network, 'Rectangle')
        p4c.export_network(filename =  artefacts_path+ "/Cytoscape/MasterRegulators_TFs_Targets_Network" , type='SIF',overwrite_file=True)
        logger.info("Network generated in Cytoscaoe software.") 
 
def create_network(artefacts_path, trrust_output,):
    logger.info(f"Creating  network using: {trrust_output}")
    merged_network_file = "MasterRegulators_TFs_Targets_Merged_Network.sif"
    generate_network(unique_node_list(trrust_output), pd.read_csv(trrust_output),artefacts_path)
    preprocess_sif(artefacts_path)
    get_regulators(get_TFs_data(artefacts_path),get_targets_data(artefacts_path),artefacts_path)
    get_master_regulators(get_TFs_data(artefacts_path),artefacts_path)
    merge_networks(artefacts_path)
    set_visual_style(merged_network_file,get_TFs_data(artefacts_path),get_targets_data(artefacts_path),artefacts_path)



    