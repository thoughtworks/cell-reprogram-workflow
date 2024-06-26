"test_functional_analysis.py"

import pytest
import shutil
import os
import filecmp
from crgem.stages.functional_analysis import functional_analysis
from crgem.common import clean_up

transsynw_genes_file = "/Trrust_Analysis/transsynw_genes.csv"
signet_genes_file = "/Trrust_Analysis/signet_genes.csv"

artefacts_path = "./test_dump"
source_artefacts_path = "./tests/artefacts"

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
    
    functional_analysis(artefacts_path,transsynw_genes_file, signet_genes_file)
    files_under_scrutiny = ['/Uniprot/uniprot_signet_genes.csv','/Uniprot/uniprot_transsynw_genes.csv' ]

    for file in files_under_scrutiny:
        file_path = file
        assert filecmp.cmp(source_artefacts_path + file_path, artefacts_path + file_path) == True
