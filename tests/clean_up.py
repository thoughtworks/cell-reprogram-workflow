"clean_up.py"
import shutil
import os

artefacts_path = "./test_dump"

def clean_up():
    shutil.rmtree(artefacts_path)
    os.mkdir(artefacts_path)