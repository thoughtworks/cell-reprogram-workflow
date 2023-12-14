"generate_hypothesis.py"


from src.tool.transsynw_tool import runTranssynW


def generate_hypothesis(artefacts_path, *args, deps_root="."):
    
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
        ARTEFACTS_PATH+"/TransSynW",
        deps_root=deps_root)
    
