"transsynw_tool.py"

import subprocess
import os

def runTranssynW(*args):
    print("RUNNING TransSynW with params-", args)
    
    STARTPOP=args[0]
    INPUT_FILE=args[1]
    INPUT_FILE2=args[2]
    SUBPOPULATIONS=args[3]
    EMAIL=args[4]
    SPECIES=args[5]
    ORGNAME=args[6]
    OUTPUT_DIR=args[7]

    os.environ["CPATH"]="/opt/homebrew/Cellar/boost/1.81.0_1/include"
    os.environ["LIBRARY_PATH"]="/opt/homebrew/lib"

    os.mkdir(OUTPUT_DIR)
    f = open(OUTPUT_DIR + "/output.txt", "w")

    script = "src/tool/code/SynergisticCore_noExpCutoff_sizeNormalization.R"
    command = "Rscript"
    command_args = [script, INPUT_FILE, SUBPOPULATIONS, OUTPUT_DIR, INPUT_FILE2, STARTPOP, SPECIES, ORGNAME]
    return subprocess.run((command, *command_args), stdout=f)



    