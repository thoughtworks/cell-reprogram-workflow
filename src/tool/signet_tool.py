"signet_tool.py"

import subprocess
import os
import sys

from src.tool.SIGNET import signet
code_path = './src/tool/Signet'

def runSignetpy(artefacts_path, ter_data):
    signet_out_path = artefacts_path + '/Signet/Intermediate'
    os.makedirs(signet_out_path)
    
    script = code_path+"/SIGNET.py"
    tf = code_path+"/tf_list.txt"
    sp = "Homo sapiens"

    old_args = sys.argv
    sys.argv = [script, "--data_file", ter_data, "--output_path", signet_out_path, "--tf_list_file", tf, "--species", sp]
    signet()
    sys.argv = old_args

def runR(artefacts_path):
    script = code_path+"/SIGNET_GRN_prediction.R"
    command = ("Rscript", script, artefacts_path)
    return subprocess.call(command)