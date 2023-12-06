"test_functional_analysis.py"

import pytest
import shutil
import os
import filecmp
from src.stages.functional_analysis import start_uniprot_transyn,start_uniprot_signet

artefacts_path = "./test_dump"
source_artefacts_path = "./tests/example"

def copy_dependencies():
    src = source_artefacts_path + '/Trrust_Analysis'
    dest = artefacts_path + '/Trrust_Analysis'
    shutil.copytree(src, dest) 

def cleanup():
    shutil.rmtree(artefacts_path, ignore_errors=True)

def setup():
    cleanup()

    os.mkdir(artefacts_path)

    copy_dependencies()

def test_should_match_generated_artefacts():
    setup()
    
    start_uniprot_transyn(artefacts_path)
    files_under_scrutiny = ['/Uniprot/transync/Uniprot_transsynw_analysis.csv' ]
    start_uniprot_signet(artefacts_path)
    files_under_scrutiny = ['/Uniprot/signet/Uniprot_signet_analysis.csv' ]

    for file in files_under_scrutiny:
        file_path = file
        assert filecmp.cmp(source_artefacts_path + file_path, artefacts_path + file_path) == True
