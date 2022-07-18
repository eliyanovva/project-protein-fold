"""This module contains the Argument Parser object with the expected arguments for modeling.
"""
import argparse

class ModelingParser(argparse.ArgumentParser):
    """This class specifies the Argument Parser which requests for the needed flags
    in modeling the ligand-protein relationship.
    """
    def __init__(self):
        super().__init__(
            description='Predict the binding affinity between your protein and ligand!',
            epilog='Thanks for stopping by!'
        )

    def setup_arguments(self):
        """Adds the expected arguments to the parser class.
        """
        #self.add_argument(
        #    '--big', help = 'Makes the side of the shape 100 pixels.', action = 'store_true'
        #)
        self.add_argument(
            '--hptuning', help = 'Try the GNN with various hyperparameters.', action = 'store_true'
        )
        self.add_argument(
            '--batch_size',
            help = 'Sets the size of the dataset to be used.',
            type = int
        )
        #self.add_argument(
        #    '--pen-color', help = 'Sets the color of the pen to a given color.'
        #)
        #self.add_argument(
        #    '--time',
        #    help = 'Sets the time to wait after the shape is completed in seconds,\
        #        defaults to 3 seconds.',
        #    type = int
        #)
        self.add_argument(
            'model',
            help = 'Specifies the type of model to be used. Choose between "cnn", "gnn", or "rf".'
        )