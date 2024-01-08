"mechanistic_insights.py"


from craft.tool.transsynw_tool import runTranssynW
from craft.tool.paga_tool import runPaga

import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
file_handler = logging.FileHandler('logs.log')
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)

def mechanistic_insights(artefacts_path, *args):

    logger.info(f"Obtaining mechanistic insights using: '{args}")

    STARTDATA=args[0]
    COMPLETEDATA=args[1]
    ANNOTATION=args[2]
    TERPOP=args[3]
    EMAIL="default@gmail.com"
    SPECIES="Human"
    TEXT=args[1]+"::"+args[2]+"::"+args[0]
    STARTPOP = args[4]
    ARTEFACTS_PATH=artefacts_path

    runTranssynW(STARTDATA, COMPLETEDATA, ANNOTATION, "::"+TERPOP+"::", EMAIL, SPECIES, TEXT, ARTEFACTS_PATH+"/TransSynW")
    runPaga(STARTPOP, TERPOP, COMPLETEDATA, ANNOTATION, ARTEFACTS_PATH)

