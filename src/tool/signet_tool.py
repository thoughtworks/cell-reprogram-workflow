"signet_tool.py"

import subprocess
import os
import sys

# from dependencies.Signet.SIGNET import signet
code_path = './Signet'

abs_path = os.path.dirname(__file__)

def runSignetpy(artefacts_path, ter_data, deps_root="."):
    signet_out_path = artefacts_path + '/Signet/Intermediate'
    os.makedirs(signet_out_path)
    # script_path = f"{abs_path}/../../deps/transsynw/"
    script = f"{abs_path}/../../deps/SIGNET.py"
    tf = code_path+"/tf_list.txt"
    sp = "Homo sapiens"

    old_args = sys.argv
    # sys.argv = [script, "--data_file", ter_data, "--output_path", signet_out_path, "--tf_list_file", tf, "--species", sp]
    command = ("python",script, "--data_file", ter_data, "--output_path", signet_out_path, "--tf_list_file", tf, "--species", sp)
    return subprocess.run(command)

    sys.argv = old_args

def runR(artefacts_path):
    script = code_path+"/SIGNET_GRN_prediction.R"
    command = ("Rscript", script, artefacts_path)
    return subprocess.call(command)