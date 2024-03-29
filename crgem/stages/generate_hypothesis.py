"generate_hypothesis.py"


from crgem.tool.transsynw_tool import runTranssynW
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
file_handler = logging.FileHandler('logs.log')
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)

def generate_hypothesis(artefacts_path, *args):
    logger.info(f"Generating hypothesis using: {args}")
    STARTDATA=args[0]
    COMPLETEDATA=args[1]
    ANNOTATION=args[2]
    TERPOP=args[3]
    EMAIL="default@gmail.com"
    SPECIES="Human"
    TEXT=args[1]+"::"+args[2]+"::"+args[0]
    ARTEFACTS_PATH=artefacts_path
    

    runTranssynW(STARTDATA,
        COMPLETEDATA,
        ANNOTATION,
        "::"+TERPOP+"::",
        EMAIL,
        SPECIES,
        TEXT,
        ARTEFACTS_PATH+"/TransSynW")
    
