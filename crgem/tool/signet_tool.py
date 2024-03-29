"signet_tool.py"

import subprocess
import os
import sys
from crgem.common import clean_up

# from dependencies.Signet.SIGNET import signet
code_path = './Signet'

abs_path = os.path.dirname(__file__)
print(abs_path)


def runSignetpy(artefacts_path, ter_data):
    clean_up(artefacts_path+'/Signet/Intermediate')
    signet_out_path = artefacts_path + '/Signet/Intermediate'
    os.makedirs(signet_out_path)
  
    script = f"{abs_path}/../dependencies/SIGNET.py"

    print("sp:",script)
    tf = f"{abs_path}/Signet/tf_list.txt"
    print("tf_path:",tf)
    sp = "Homo sapiens"

    old_args = sys.argv
    command = ("python",script, "--data_file", ter_data, "--output_path", signet_out_path, "--tf_list_file", tf, "--species", sp)
    print(command)
    return subprocess.run(command)


def runR(artefacts_path):
    script =f"{abs_path}/Signet/SIGNET_GRN_prediction.R"
    command = ("Rscript", script, artefacts_path)
    return subprocess.call(command)