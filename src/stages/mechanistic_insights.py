"mechanistic_insights.py"

from src.tool.transsynw_tool import runTranssynW
from src.tool.paga_tool import runPaga

def mechanistic_insights(artefacts_path, *args):
    print("RUNNING stage 'mechanistic_insights' with params-", args)

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

