#This script imports TM domain sequences and positions from the Uniprot database

import json

#Create comma-separated list of accession numbers
with open('all_accessions.txt') as f: #TODO: Try with all accessions
    lines = f.readlines()

query = []
i = 0
addit = ""
for line in lines:
  if i < 100:
    addit += line.replace('\n','') + "%2C"
  else:
    query.append(addit)
    i = 0
    addit = ""
  if (line == lines[-1]) & i!=0:
    query.append(addit)
  i+=1
  

#The following code uses the uniprot API: https://www.ebi.ac.uk/proteins/api/doc/#/ to lookup transmembrane features for all proteins
#https://academic.oup.com/nar/article/45/W1/W539/3106040?login=false

import requests, sys

with open("TM.txt", "w") as f:
  for set in query:
    requestURL = "https://www.ebi.ac.uk/proteins/api/features?offset=0&size=100&accession=" + set + "&types=TRANSMEM"

    r = requests.get(requestURL, headers={ "Accept" : "application/json"})

    if not r.ok:
      r.raise_for_status()
      sys.exit()

    #Results from scraping
    responseBody = r.text
    TM_dict = json.loads(responseBody)
    #print(json.dumps(TM_dict, indent=4))
    for accession in TM_dict:
      if len(accession["features"]) > 6:
        print(accession['accession'], file=f)
        print(accession['sequence'][int(accession["features"][2]['begin'])-1:int(accession["features"][2]['end'])], file=f)
        print(accession['sequence'][int(accession["features"][4]['begin'])-1:int(accession["features"][4]['end'])], file=f)
        print(accession['sequence'][int(accession["features"][5]['begin'])-1:int(accession["features"][5]['end'])], file=f)
        print(accession['sequence'][int(accession["features"][6]['begin'])-1:int(accession["features"][6]['end'])], file=f)
