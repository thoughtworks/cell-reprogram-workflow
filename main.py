
import logging
import sys
import os
from src.stages.functional_analysis import start_uniprot_transyn, start_uniprot_signet
from src.stages.mechanistic_insights import mechanistic_insights
from src.stages.grn_inference import grn_inference
from src.stages.trrust_analysis import trrust_analysis
from src.stages.create_network import create_network
from src.stages.generate_hypothesis import generate_hypothesis

log_file = open("terminal_output.log", "w")
# sys.stdout = log_file
# sys.stderr = log_file

# Configure the basic logging
logging.basicConfig(level=logging.INFO)

# logger = logging.getLogger()
# logger.setLevel(logging.INFO)
# formatter = logging.Formatter('%(asctime)s | %(levelname)s | %(message)s')
# file_handler = logging.FileHandler('logs.log')
# file_handler.setLevel(logging.DEBUG)
# file_handler.setFormatter(formatter)
# logger.addHandler(file_handler)

# logger.info("usage: run stage-name --artefacts <artefacts_path> --params [...params]")
# logger.info("available stages: transsynw, mechanistic_insights, grn_inference, trrust_analysis, generate_network, functional_analysis")
       
abs_path = os.path.dirname(__file__)
DEPS_PATH = f"{abs_path}/deps"


def run_all(artefacts_path, *args):
    TERDATA = args[5]
    TRRUST_DB_FILE = args[6] or "./trrust_rawdata_human.tsv"
    

    mechanistic_insights(artefacts_path, *args)
    grn_inference(artefacts_path=artefacts_path, ter_data=TERDATA)
    trrust_analysis(TRRUST_DB_FILE, artefacts_path)
    create_network(artefacts_path + "/Trrust_Analysis/trrust_analysis.csv",artefacts_path)
    start_uniprot_transyn(artefacts_path)
    start_uniprot_signet(artefacts_path)


def main():
    stage_name = sys.argv[2]
    artefacts_path = os.path.abspath(sys.argv[4])
    # convert string to absolute path
    
    if stage_name == "generate_hypothesis":
        params = sys.argv[6:]
        generate_hypothesis(artefacts_path, *params, deps_root=DEPS_PATH)
    elif stage_name == "mechanistic_insights":
        params = sys.argv[6:]
        mechanistic_insights(artefacts_path, *params)
    elif stage_name == "grn_inference":
        TERDATA = sys.argv[6]
        grn_inference(deps_root=DEPS_PATH, artefacts_path=artefacts_path, ter_data=TERDATA)
    elif stage_name == "trrust_analysis":
        trrust_db_file = sys.argv[6]
        trrust_analysis(trrust_db_file, artefacts_path)
    elif stage_name == "create_network":
        
        trrust_output = sys.argv[6]
        create_network(trrust_output, artefacts_path)
    elif stage_name == "functional_analysis":
        start_uniprot_transyn(artefacts_path)
        start_uniprot_signet(artefacts_path)
    elif stage_name == "all":
        run_all(artefacts_path, *sys.argv[6:])
    else:
        print("usage: run stage-name --artefacts <artefacts_path> --params [...params]")
        print("available stages: transsynw, mechanistic_insights, grn_inference, trrust_analysis, generate_network, functional_analysis")
        sys.exit(1)

main()