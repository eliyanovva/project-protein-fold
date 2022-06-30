import numpy as np
import math
import os
import constants
#import logging as log
#import log_config
from feature_matrices import getMatrix

class PDBDataFile:
    def __init__(self, pdb_filename):
        #log.info('class initialized :)')
        self.pdb_filename = pdb_filename
    
    def getFeatureMatrix(self):
        getMatrix(self.pdb_filename)