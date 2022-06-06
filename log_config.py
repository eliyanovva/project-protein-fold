import logging
import sys

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
