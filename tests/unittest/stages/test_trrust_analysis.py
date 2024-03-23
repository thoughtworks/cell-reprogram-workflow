"test_trrust_analysis.py"

import pytest
import shutil
import os
import filecmp
from cmi.stages.trrust_analysis import trrust_analysis

trrust_db_file = "./data/trrust_rawdata_human.tsv"
artefacts_path = "./test_dump"
source_artefacts_path = "./tests/artefacts"

def copy_dependencies():
    src = source_artefacts_path + '/TransSynW'
    dest = artefacts_path + '/TransSynW'
    shutil.copytree(src, dest) 

    src = source_artefacts_path + '/Signet'
    dest = artefacts_path + '/Signet'
    shutil.copytree(src, dest) 

def cleanup():
    shutil.rmtree(artefacts_path, ignore_errors=True)
    
def setup():
    cleanup()

    os.makedirs(artefacts_path)
    copy_dependencies()

def test_should_match_generated_artefacts():
    setup()

    trrust_analysis(trrust_db_file, artefacts_path)
    files_under_scrutiny = ["trrust_analysis.csv", "transsynw_genes.csv", "signet_genes.csv" ]

    for file in files_under_scrutiny:
        file_path =  "/Trrust_Analysis/" + file
        assert filecmp.cmp(source_artefacts_path + file_path, artefacts_path + file_path) == True
    