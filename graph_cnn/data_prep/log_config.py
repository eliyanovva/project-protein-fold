import logging
import sys
import os

# Create logging directory
if not os.path.exists('./logs'):
    os.makedirs('./logs')

    
logging.basicConfig(
    encoding='utf-8', 
    level=logging.INFO, 
    format='%(asctime)s %(message)s',
    filename='logs/script.log',

)
