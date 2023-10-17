"grn_inference.py"

from src.tool.signet_tool import runSignetpy, runR


def grn_inference(**kwargs):
    print("RUNNING grn_inference with params", kwargs)

    TERDATA=kwargs['ter_data']
    ARTEFACTS_PATH = kwargs['artefacts_path']

    runSignetpy(ARTEFACTS_PATH, TERDATA)
    runR(ARTEFACTS_PATH)