"transsynw_tool.py"

import logging
import subprocess
import os
import sys
from craft.common import clean_up
cwd = os.getcwd()
print("cwd:", cwd)
abs_path = os.path.dirname(__file__)
print("abs_path", abs_path)

def runTranssynW(*args):
    STARTPOP=args[0]
    INPUT_FILE=args[1]
    INPUT_FILE2=args[2]
    SUBPOPULATIONS=args[3]
    EMAIL=args[4]
    SPECIES=args[5]
    ORGNAME=args[6]
    OUTPUT_DIR=args[7]
    

    clean_up(OUTPUT_DIR)
    os.mkdir(OUTPUT_DIR)
    f = open(OUTPUT_DIR + "/output.txt", "w")
    present_dir = os.getcwd()
    data_dir_path = f"{cwd}/data/"
    print("data_path", data_dir_path)
    script_path = f"{abs_path}/../dependencies/transsynw/"
    print("script_path", script_path)
    os.chdir(script_path)
    script = "code/SynergisticCore_noExpCutoff_sizeNormalization.R"
    command = "Rscript"
    command_args = [script, data_dir_path + INPUT_FILE, SUBPOPULATIONS, OUTPUT_DIR, data_dir_path + INPUT_FILE2, data_dir_path + STARTPOP, SPECIES, ORGNAME]
    return_val = subprocess.run((command, *command_args), stdout=f)
    os.chdir(present_dir)

    return return_val






    