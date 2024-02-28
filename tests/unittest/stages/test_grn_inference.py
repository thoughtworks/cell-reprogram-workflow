"test_grn_inference.py"

import pytest
import shutil
import os
import filecmp
from craft.stages.grn_inference import runR
# from craft.common import clean_up

artefacts_path = "./test_dump"
source_artefacts_path = "./tests/example"

def copy_dependencies():

    src = source_artefacts_path + '/Signet/Intermediate'
    dest = artefacts_path + '/Signet/Intermediate'
    shutil.copytree(src, dest) 

def cleanup():
    shutil.rmtree(artefacts_path, ignore_errors=True)

def setup():
    cleanup()
    os.mkdir(artefacts_path)
    copy_dependencies()

def test_should_match_generated_artefacts():
    setup()
    
    runR(artefacts_path)
    files_under_scrutiny = ['/Signet/copaired2.csv']

    for file in files_under_scrutiny:
        file_path = file
        assert filecmp.cmp(source_artefacts_path + file_path, artefacts_path + file_path) == True
