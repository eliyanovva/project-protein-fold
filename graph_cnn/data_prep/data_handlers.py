from abc import abstractmethod
import numpy as np
import pandas as pd
import tensorflow as tf
import logging as log
import os
import sys
import re

from . import log_config
from . import constants
from abc import ABC, abstractmethod

class DataHandlers(ABC):
    # THE CLASS IS MADE TO WORK ONLY WITH 2D MATRICES
    def __init__(self):
        pass


    def initialize(self, labels_list, dimensions, type):
        # let dimensions be a tuple, not a list, so it can be passed as an argument
        # to numpy and tensorflow functions.
        self.labels_list = labels_list
        self.dimensions = dimensions
        self.type = type
        log.info('Initialized data handler class')
        

    def __fillMatrix(self):
        log.info('Started filling out data matrix')
        for i in range(len(self.labels_list)):
            data = self.loadDataSingleMatrix(self.labels_list[i])
            log.info('Succesfully loaded data for ' + self.labels_list[i])
            self.matrix[i] = np.pad(
                np.array(data, ndmin=2),
                ((0, self.dimensions[0] - len(data)),
                (0, self.dimensions[1] - len(data[0]))),
                'constant',
                constant_values=(0)
            )
            log.info('Padded and added data for ' + self.labels_list[i])

    
    def getTensor(self):
        log.info('Started Tensor extraction process')
        self.__initializeNpMatrix()
        self.__fillMatrix()
        return tf.convert_to_tensor(self.matrix)

    @abstractmethod
    def loadDataSingleMatrix(self, label_name):
        raise NotImplementedError("Must override methodB")

    
    def __initializeNpMatrix(self):
        self.matrix = np.zeros(
            (
                len(self.labels_list),
                self.dimensions[0],
                self.dimensions[1]
            ),
            dtype = self.type
        )
        log.info('Completed zero matrix initialization')
    