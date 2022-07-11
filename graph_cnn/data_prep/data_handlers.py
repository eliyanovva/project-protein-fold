"""The DataHandlers class helps generate numpy arrays which are fed into the GNN model.
It uses the .npy adjacency and feature matrices which are already generated. Its unique
features are that it makes them the same size, and generated a big numpy array for each
type of matrix, which can hold all training protein adjacency data, ot all training feature data
and so on.
"""
import logging as log
from abc import ABC, abstractmethod

import numpy as np
import tensorflow as tf

import config

class DataHandlers(ABC):
    """This is an abstract class used to generate Tensors from .npy matrices.
    It is made to work only with 2D matrices.

    Raises:
        NotImplementedError: This error is raised if the method used to extract
        specific matrix type from a .npy file is not defined.
    """
    def __init__(self):
        pass


    def initialize(self, labels_list, dimensions, data_type):
        """Initialization method for the DataHandlers class.

        Args:
            labels_list (List[str]): A list of strings which are protein or ligand names.
            dimensions (tuple): A tuple (x, y) with the dimensions of the matrix type to be handled
            type (np.type): the datatype of the matrix - int or float.
        """
        self.labels_list = labels_list
        self.dimensions = dimensions
        self.type = data_type
        log.info('Initialized data handler class')


    def getTensor(self):
        """Returns a tensor containing all requested matrices.

        Returns:
            tf.Tensor: A Tensor which concatenates all matrices labeled in the labels_list.
        """
        log.info('Started Tensor extraction process')
        self.__initializeNpMatrix()
        self.__fillMatrix()
        return tf.convert_to_tensor(self.matrix)

    @abstractmethod
    def loadDataSingleMatrix(self, label_name):
        """This method loads a single .npy array to a numpy array and returns it. It has
        to be defined be every child class.

        Args:
            label_name (string): The name of the protein/ligand whose data we need.

        Raises:
            NotImplementedError: Raised if the method is not implemented in the child class.
        """
        raise NotImplementedError("Must override methodB")


    def __fillMatrix(self):
        log.info('Started filling out data matrix')
        for i in range(len(self.labels_list)):
            data = self.loadDataSingleMatrix(self.labels_list[i])
            log.info('Succesfully loaded data for ' + self.labels_list[i])
            data_rows = len(data)
            data_cols = len(data[0])
            
            if data_rows < self.dimensions[0] and data_cols < self.dimensions[1]:
                self.matrix[i] = np.pad(
                    np.array(data, ndmin=2),
                    ((0, self.dimensions[0] - data_rows),
                    (0, self.dimensions[1] - data_cols)),
                    'constant',
                    constant_values=(0)
                )
                log.info('Padded and added data for ' + self.labels_list[i])
          
            elif data_rows > self.dimensions[0] or data_cols > self.dimensions[1]:
                row_diff = (data_rows - self.dimensions[0]) // 2
                col_diff = (data_cols - self.dimensions[1]) // 2
                #print(data_rows - self.dimensions[0])
                self.matrix[i] = np.array(data, ndmin=2)[
                    row_diff : self.dimensions[0] + row_diff ,
                    col_diff : self.dimensions[1] + col_diff 
                ]
                log.info('Cropped and added data for ' + self.labels_list[i])

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
