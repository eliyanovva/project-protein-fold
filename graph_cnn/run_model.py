from lib2to3.pytree import convert
import time
from graph_cnn.model import GraphCNN
import config
import logging as log
import numpy as np
import tensorflow as tf
#import visuals

def convertYClassification(y):
    new_y = np.zeros(len(y), dtype=bool)
    for i in range(len(y)):
        new_y[i] = abs(y[i]) >= 1
    return tf.convert_to_tensor(new_y, dtype=bool)


def createCallbacks():
    callbacks = [
        tf.keras.callbacks.EarlyStopping(
            monitor = 'loss',#was val_loss before
            min_delta=0.05,
            patience = 0,
            restore_best_weights = True,
            mode = 'min',
        ),
        tf.keras.callbacks.TensorBoard(
            log_dir='logs/tensorboard',
            histogram_freq=0,
            write_graph=True,
            write_images=False,
            write_steps_per_second=False,
            update_freq='epoch',
            profile_batch=0,
            embeddings_freq=0,
            embeddings_metadata=None,
        ),
    ]

    return callbacks


def trainGNN(gnn, X_train, y_train, tr_batch_size=32, tr_optimizer='adagrad', classification=False):
    if classification:
        model = gnn.classificationModel(hp_optimizer=tr_optimizer)
        y_train = convertYClassification(y_train)
    else:
        model = gnn.createModel(hp_optimizer=tr_optimizer)
    callbacks = createCallbacks()
    log.info('model fitting started')
    mod_history = model.fit(
        X_train, y_train, epochs=10, verbose=True,
        batch_size=tr_batch_size, callbacks = callbacks, validation_split=0.2)
    
    log.info ('model fitting finished successfully')

    return model, mod_history


def testGNN(model, X_test, y_test, classification=False):
    log.info('model evaluation started')
    if classification:
        y_test = convertYClassification(y_test)
    results = model.evaluate(X_test, y_test, verbose=True)
    log.info('The results of the test run are:' + ' '.join([str(r) for r in results]))
    log.info('model evaluation completed')
    return results


def runGNN(model, X_run):
    """Runs the trained model with new data.

    Args:
        model (tf.keras.Model): The pretrained GNN model.
        X_run (List[tf.Tensor]): The 4 tensors corresponding to the input data
    """
    return model.predict(X_run)

def storeResults(results, timing_measures):
    with open('results.txt', 'a') as res_log:
        res_log.write(' '.join([str(r) for r in results]) + ' \n')
        res_log.write('Timing Benchmarks:\n')
        res_log.write(' '.join([str(r) for r in timing_measures]) + '\n')


def evaluateTuple(model, X_val):
    """Should predict the binding coefficient of a new tuple.

    Args:
        model (_type_): _description_
        X_val (_type_): _description_
    """
    pass


def runModel(batch_size=-1, test_frac=0.3, classification=False):
    gnn = GraphCNN()
    gnn.initialize()
    
    start_train_test_split = time.time()
    X_train_labels, X_test_labels, y_train_labels, y_test_labels = gnn.trainTestSplit(
        model_test_size=test_frac, batch_size=batch_size)
    end_train_test_split = time.time()

    start_train_data_load = time.time()
    X_train, y_train = gnn.getTensors(X_train_labels, y_train_labels)
    end_train_data_load = time.time()

    start_test_data_load = time.time()
    X_test, y_test = gnn.getTensors(X_test_labels, y_test_labels)
    end_test_data_load = time.time()

    start_model_fitting = time.time()
    model, training_history = trainGNN(gnn, X_train, y_train, classification=classification)
    end_model_fitting = time.time()

    start_model_fitting = time.time()
    results = testGNN(model, X_test, y_test, classification=classification)
    end_model_fitting = time.time()

    timing_measures = [
        end_train_test_split - start_train_test_split,
        end_train_data_load - start_train_data_load,
        end_test_data_load - start_test_data_load,
        end_model_fitting - start_model_fitting,
        end_model_fitting - start_train_test_split
    ]

    storeResults(results, timing_measures)

    
    return model


def evalNewData(new_X):
    pass
#runModel()
