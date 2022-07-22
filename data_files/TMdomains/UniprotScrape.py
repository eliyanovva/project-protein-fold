#This script imports TM domain sequences and positions from the Uniprot database

import json

def scrape_TMs(proteins, writefile, csv):
  #Create comma-separated list of accession numbers
  with open(proteins) as f: #TODO: Try with all accessions
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

  with open(writefile, "w") as f:
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
          print(('>' + accession['accession']), file=f)
          print((accession['sequence'][int(accession["features"][2]['begin'])-6:int(accession["features"][2]['end'])+5]+','+str(int(accession["features"][2]['begin'])-5)+','+str(int(accession["features"][2]['end'])+5)), file=f)
          print((accession['sequence'][int(accession["features"][4]['begin'])-6:int(accession["features"][4]['end'])+5]+','+str(int(accession["features"][4]['begin'])-5)+','+str(int(accession["features"][4]['end'])+5)), file=f)
          print((accession['sequence'][int(accession["features"][5]['begin'])-6:int(accession["features"][5]['end'])+5]+','+str(int(accession["features"][5]['begin'])-5)+','+str(int(accession["features"][5]['end'])+5)), file=f)
          print((accession['sequence'][int(accession["features"][6]['begin'])-6:int(accession["features"][6]['end'])+5]+','+str(int(accession["features"][6]['begin'])-5)+','+str(int(accession["features"][6]['end'])+5)), file=f)
  
  #Convert to csv
  with open(writefile) as f:
    lines = f.readlines()

  with open(csv, 'w') as f:
      printstatement = 'protein,TM3,s3,e3,TM5,s5,e5,TM6,s6,e6,TM7,s7,e7'
      for line in lines:
          if line[0] == '>':
              print('stop')
              print(printstatement)
              print(printstatement, file = f)
              print('start')
              line = line.replace('>', '')
              printstatement = line.replace('\n', '')
              print(printstatement)
          else:
              printstatement += ',' + line.replace('\n', '')

      print(printstatement, file = f)

scrape_TMs('all_accessions.txt', 'TM.txt', 'TM.csv')