#This script imports TM domain sequences and positions from the Uniprot database

import json

#Create comma-separated list of accession numbers
with open('accessions_tolookup.txt') as f: #TODO: Try with all accessions
    lines = f.readlines()

query = ""
for line in lines:
    query += line.replace('\n','')
    query+= "%2C"

#The following code uses the uniprot API: https://www.ebi.ac.uk/proteins/api/doc/#/ to lookup transmembrane features for all proteins
#https://academic.oup.com/nar/article/45/W1/W539/3106040?login=false

import requests, sys

requestURL = "https://www.ebi.ac.uk/proteins/api/features?offset=0&size=100&accession=" + query + "&types=TRANSMEM"

r = requests.get(requestURL, headers={ "Accept" : "application/json"})

if not r.ok:
  r.raise_for_status()
  sys.exit()

#Prints results from scraping
responseBody = r.text
TM_dict = json.loads(responseBody)
#print(json.dumps(TM_dict, indent=4))

for accession in TM_dict:
  print(accession['accession'])
