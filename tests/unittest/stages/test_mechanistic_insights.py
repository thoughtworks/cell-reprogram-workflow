"test_mechanistic_insights.py"

import filecmp
import os
import shutil
import glob
from crgem.stages.mechanistic_insights import mechanistic_insights

args = ["start.txt","start_terminal.txt","annotation.txt","HPROGFPM", "HNES"]
artefacts_path = os.getcwd() + "/" +  "./test_dump"
source_artefacts_path = "./tests/artefacts"

def cleanup():
    shutil.rmtree(artefacts_path)

def setup():
    cleanup()

    os.mkdir(artefacts_path)

def test_should_match_generated_artefacts():
    setup()
    
    mechanistic_insights(artefacts_path, *args)
    files_under_scrutiny = ["/TransSynW/cores.tsv", "/TransSynW/markers.tsv"]
    images = glob.glob("./test_dump/PAGA/*.png", recursive=True)


    for file in files_under_scrutiny:
        # file_path =  "/TransSynW/" + file
        assert filecmp.cmp(source_artefacts_path + file, artefacts_path + file) == True
        
    assert len(images) != 0
    

