"grn_inference.py"

from src.tool.signet_tool import runSignetpy, runR


def grn_inference(deps_root=".", **kwargs):
    print("RUNNING grn_inference with params", kwargs)

    TERDATA=kwargs['ter_data']
    ARTEFACTS_PATH = kwargs['artefacts_path']

    runSignetpy(ARTEFACTS_PATH, TERDATA, deps_root=deps_root)
    runR(ARTEFACTS_PATH)