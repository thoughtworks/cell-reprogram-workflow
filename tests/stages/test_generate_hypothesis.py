"test_generate_hypothesis.py"

import filecmp
import os
import shutil
from src.stages.generate_hypothesis import generate_hypothesis

args = ["start.txt","start_terminal.txt","annotation.txt","HPROGFPM", "HNES"]
artefacts_path = "./test_dump"
source_artefacts_path = "./tests/artefacts"

print("Executed")

def cleanup():
    shutil.rmtree(artefacts_path)

def setup():
    cleanup()

    os.mkdir(artefacts_path)




def test_should_match_generated_artefacts():
    setup()

    generate_hypothesis(artefacts_path, *args)
    files_under_scrutiny = ["cores.tsv", "markers.tsv"]

    for file in files_under_scrutiny:
        file_path =  "/TransSynW/" + file
        assert filecmp.cmp(source_artefacts_path + file_path, artefacts_path + file_path) == True
    



        

