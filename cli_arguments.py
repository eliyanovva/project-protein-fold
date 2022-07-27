"""This module contains the Argument Parser object with the expected arguments for modeling.
"""
import argparse
import string

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
        self.add_argument(
            '--gnn_cl', help = 'Runs the GNN as a binary classificator instead of a regressor.', action = 'store_true'
        )

        self.add_argument(
            '--batch_size',
            help = 'Sets the size of the dataset to be used.',
            type = int
        )

        self.add_argument(
            '--optimizer',
            help = 'Sets the optimizer.',
            type = string
        )

        self.add_argument(
            '--dropout',
            help = 'Sets the size of the dropout.',
            type = float
        )

        self.add_argument(
            '--test_train_split',
            help = 'Sets the size of the test train split.',
            type = float
        )

        self.add_argument(
            '--validation_split',
            help = 'Sets the validation split.',
            type = float
        )

        self.add_argument(
            '--learning_rate',
            help = 'Sets the learning rate. Adjusts depending on optimizer.',
            type = float
        )

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

        # run - splits the set into training/testing, and runs the model on that
        # hptuning - runs the hptuning script
        # eval_tuple - takes a protein and ligand, converts them to matrices and runs them through the model
        # eval_protein - takes a protein and evaluates its binding coefficients with all available ligands
        # eval_ligand - takes a ligand and evaluates its binding coefficient with all available proteins.
        self.add_argument(
            '--gnn_mode',
            help = 'Choose between "run", "hptuning", "eval_tuple", "eval_ligand", "eval_protein". The program defaults to the "run" mode.'
        )

        self.add_argument(
            '--interaction',
            help = "Input interaction coefficient between the protein and ligand inputs."
        )

        self.add_argument(
            '--verbose',
            help = "If used, the program will output and preserve DEBUG level logs."
        )

        self.add_argument(
            '--rf_mode',
            help = 'Choose between "run", "eval_pairs", "eval_ligands", "eval_proteins". The program defaults to the "run" mode.'
        )        