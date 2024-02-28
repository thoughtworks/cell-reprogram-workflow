import filecmp
import pytest
import os
import shutil
import glob
import pandas as pd
from craft.stages.create_network import create_network #, get_TFs_data, get_targets_data, get_regulators, get_master_regulators


artefacts_path = "./test_dump"
source_artefacts_path = "./tests/artefacts"
trrust_output = artefacts_path + '/Trrust_analysis/trrust_analysis.csv'

def copy_dependencies():
    src = source_artefacts_path + '/Trrust_analysis'
    dest = artefacts_path + '/Trrust_analysis'
    shutil.copytree(src, dest) 

    src = source_artefacts_path + '/TransSynW'
    dest = artefacts_path + '/TransSynW'
    shutil.copytree(src, dest) 

    src = source_artefacts_path + '/Signet'
    dest = artefacts_path + '/Signet'
    shutil.copytree(src, dest) 


def cleanup():
    shutil.rmtree(artefacts_path)

def setup():
    cleanup()

    os.makedirs(artefacts_path)
    copy_dependencies()
    
def test_TFs():
    setup()
    create_network(artefacts_path, trrust_output)
    files_under_scrutiny = ["MasterRegulators_TFs_Targets_Merged_Network.sif", "PathLinker_MasterRegulators_TFs_Paths.sif", "PathLinker_TF_Targets_Paths.sif" ]

    for file in files_under_scrutiny:
        file_path =  "/Cytoscape/" + file
        assert filecmp.cmp(source_artefacts_path + file_path, artefacts_path + file_path) == True

    