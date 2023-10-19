
import sys
from src.stages.functional_analysis import start_uniprot
#from src.common import unique_directory_name
from src.stages.mechanistic_insights import mechanistic_insights
from src.stages.grn_inference import grn_inference
from src.stages.trrust_analysis import trrust_analysis, combine_transync_outputs
from src.stages.generate_network import generate_network
from src.stages.generate_hypothesis import generate_hypothesis

def run_all(artefacts_path, *args):
    TERDATA = args[5]
    TRRUST_DB_FILE = args[6] or "./trrust_rawdata_human.tsv"

    mechanistic_insights(artefacts_path, *args)
    grn_inference(artefacts_path=artefacts_path, ter_data=TERDATA)
    trrust_analysis(TRRUST_DB_FILE, artefacts_path)
    generate_network(artefacts_path + "/trrust_analysis.csv")
    start_uniprot(artefacts_path)


def main():
    stage_name = sys.argv[2]
    artefacts_path = sys.argv[4]
    
    if stage_name == "generate_hypothesis":
        params = sys.argv[6:]
        generate_hypothesis(artefacts_path, *params)
    elif stage_name == "mechanistic_insights":
        params = sys.argv[6:]
        mechanistic_insights(artefacts_path, *params)
    elif stage_name == "grn_inference":
        TERDATA = sys.argv[6]
        grn_inference(artefacts_path=artefacts_path, ter_data=TERDATA)
    elif stage_name == "trrust_analysis":
        trrust_db_file = sys.argv[6]
        trrust_analysis(trrust_db_file, artefacts_path)
    elif stage_name == "generate_network":
         trrust_analysis_out = sys.argv[7]
         generate_network(trrust_analysis_out)
    elif stage_name == "functional_analysis":
         start_uniprot(artefacts_path)
    elif stage_name == "all":
        run_all(artefacts_path, *sys.argv[6:])
    else:
        print("usage: run stage-name --artefacts <artefacts_path> --params [...params]")
        print("available stages: transsynw, mechanistic_insights, grn_inference, trrust_analysis, generate_network, functional_analysis")
        sys.exit(1)

main()