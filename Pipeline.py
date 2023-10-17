
import os
import sys
import subprocess
import pandas as pd
import py4cytoscape as p4c
import datetime
import requests
import json

mydir = os.path.join(os.getcwd(), "Results-"+datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
os.makedirs(mydir)

STARTDATA=sys.argv[1]
COMPLETEDATA=sys.argv[2]
ANNOTATION=sys.argv[3]
TERPOP=sys.argv[4]
EMAIL="default@gmail.com"
SPECIES="Human"
TEXT=sys.argv[2]+"::"+sys.argv[3]+"::"+sys.argv[1]
STARTPOP = sys.argv[5]
TERDATA = sys.argv[6]
OUTPUTDIR = mydir

# independent
def runTranssyn(*args):
    script = "./transsyn_wrapper.sh"
    command = (script, *args)
    return subprocess.run(command)

runTranssyn(STARTDATA, COMPLETEDATA, ANNOTATION, "::"+TERPOP+"::", EMAIL, SPECIES, TEXT, OUTPUTDIR+"/TransSynW")

def runPython(*args):
    script = "paga_script.py"
    command = ("python",script, *args)
    return subprocess.run(command)

runPython(STARTPOP, TERPOP, COMPLETEDATA, ANNOTATION, OUTPUTDIR)


sigpath = mydir+'/Signet/Intermediate'
os.makedirs(sigpath)

def runSignetpy():
    script = "./Signet/SIGNET.py"
    outpath = mydir+"/Signet/Intermediate"
    tf = "./Signet/tf.txt"
    sp = "Homo sapiens"
    command = ("python",script, "--data_file", TERDATA, "--output_path", outpath, "--tf_list_file", tf, "--species", sp)
    return subprocess.run(command)

runSignetpy()


def runR():
    command = ("Rscript", "./Signet/final_signet.R", OUTPUTDIR)
    return subprocess.call(command)

runR()

#subprocess.call("Rscript ./Signet/final_signet.R", shell=True)

#TRRUST analysis
trrust_data = pd.read_table("./trrust_rawdata_human.tsv",header=None)
trrust_data.columns =['Gene', 'Target', 'Action', 'Refrence']
trrust_data = trrust_data.drop(['Refrence'],axis = 1)


transync = pd.read_table(mydir+"/TransSynW/cores.tsv")
transync.rename(columns = {'core':'Type'}, inplace = True)
transynm = pd.read_table(mydir+"/TransSynW/markers.tsv")
transynm.insert(3, 'Type', 'marker')
transyn_data = pd.concat([transync, transynm])
#transyn_data2=transyn_data.drop(['FoldChange','Subpopulation','scJSD'],axis = 1)

signet = pd.read_csv(mydir+"/Signet/copaired2.csv")
tf = signet['V1']
tg = signet['V2']
signet_data = pd.concat([tf,tg])
signet_data=pd.DataFrame(signet_data,columns=['Gene'])

ttg=trrust_data[trrust_data['Gene'].isin(transyn_data['Gene'])]
ttg['Source']='TranSyn'
ttt=trrust_data[trrust_data['Target'].isin(transyn_data['Gene'])]
ttt['Source']='TranSyn'
tsg=trrust_data[trrust_data['Gene'].isin(signet_data['Gene'])]
tsg['Source']='SIGNET'
tst=trrust_data[trrust_data['Target'].isin(signet_data['Gene'])]
tst['Source']='SIGNET'
trrust_analysis = [ttg,ttt,tsg,tst]
tdata= pd.concat(trrust_analysis)
tdata.to_csv("./trrust_analysis.csv",index=False)
a = pd.read_csv("./trrust_analysis.csv")

G = a['Gene']
T = a['Target']
merge = [G,T]
merge = pd.concat(merge)
merge.to_list()
m = merge.unique()


# Cytoscape

p4c.cytoscape_ping()
p4c.cytoscape_version_info()
tnodes = pd.DataFrame(data={'id': m})
print(tnodes)
tedges = pd.DataFrame(data={'source': tdata['Gene'], 'target': tdata['Target'], 'interaction': tdata['Action']})
print(tedges)
p4c.create_network_from_data_frames(tnodes, tedges, title="Network", collection="Trrust_Database")
p4c.set_visual_style('Marquee')


#Uniprot
ts=[]
for i in transyn_data['Gene']:
    ts.append(i)

WEBSITE_API = "https://rest.uniprot.org/"
def get_url(url, **kwargs):
    response = requests.get(url, **kwargs);
    if not response.ok:
        print(response.text)
        response.raise_for_status()
        sys.exit()
    return response
jobid=[]
for i in range(len(ts)):
    r = requests.post(f"{WEBSITE_API}/idmapping/run", data={"from": "Gene_Name","to": "UniProtKB","ids": ts[i]})
    jobid.append(r.json()['jobId'])

uni=[]
for j in jobid:
    r = get_url(f"{WEBSITE_API}/idmapping/status/{j}")
    data = r.json()
    for results in data['results']:
        try:
            if results['to']['organism']['scientificName'] == "Homo sapiens":
                uni.append(([results['from'], results['to']['entryType'], results['to']['primaryAccession'], results['to']['comments'][0]['texts'][0]['value']]))
        except:
          print("An exception occurred")
uni = pd.DataFrame(uni, columns=['Gene','Entry_type','Primary_acc','Function'])
os.mkdir(OUTPUTDIR+"/Uniprot")
uni.to_csv(OUTPUTDIR+"/Uniprot/Uniprot_ts_analysis.csv", index= False)
