import logging
import sys
import os

# Create logging directory
if not os.path.exists('./logs'):
    os.makedirs('./logs')

    
logging.basicConfig(
    #filename='logs/script.log', 
    encoding='utf-8', 
    level=logging.INFO, 
    format='%(asctime)s %(message)s',
    handlers=[
        logging.FileHandler("logs/script.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
