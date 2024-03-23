"grn_inference.py"

from cmi.tool.signet_tool import runSignetpy, runR
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
file_handler = logging.FileHandler('logs.log')
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)


def grn_inference(deps_root=".", **kwargs):
    logger.info("Inferring Gene regulatory network information using: ", kwargs)

    TERDATA=kwargs['ter_data']
    ARTEFACTS_PATH = kwargs['artefacts_path']

    runSignetpy(ARTEFACTS_PATH, TERDATA)
    runR(ARTEFACTS_PATH)