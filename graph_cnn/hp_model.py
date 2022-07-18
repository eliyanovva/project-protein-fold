#import os
#import sys

# FIXME: VERY UGLY FIX FOR THE config.py
# If the file is called through entry_points, this can be removed.
#MAIN_PACKAGE_DIR = os.path.abspath(os.curdir)
#last_dir = os.path.split(os.curdir)[-1]
#while last_dir != "project-protein-fold":
#    os.chdir("..")
#    MAIN_PACKAGE_DIR = os.path.abspath(os.curdir)
#    last_dir = os.path.split(MAIN_PACKAGE_DIR)[-1]
    
#sys.path.append(MAIN_PACKAGE_DIR)

from graph_cnn.model import GraphCNN
import graph_cnn.run_model as run_model
import time
import tensorflow as tf
import config
from tensorboard.plugins.hparams import api as hp
import logging as log


with tf.summary.create_file_writer('logs/hparam_tuning').as_default():
    hp.hparams_config(
        hparams=[config.HP_BATCH_SIZE, config.HP_DROPOUT, config.HP_OPTIMIZER],
        metrics=[hp.Metric(config.METRIC_ACCURACY, display_name='Accuracy')],
    )

def hpBuildModel(hparams = {
                #config.HP_NUM_UNITS: num_units,
                #config.HP_DROPOUT: dropout_rate,
                config.HP_BATCH_SIZE: 32,
                config.HP_OPTIMIZER: 'adagrad',
            }):
    g = GraphCNN()
    g.initialize()

    start_train_test_split = time.time()
    X_train_labels, X_test_labels, y_train_labels, y_test_labels = g.trainTestSplit()
    end_train_test_split = time.time()

    start_train_data_load = time.time()
    X_train, y_train = g.getTensors(X_train_labels, y_train_labels)
    end_train_data_load = time.time()

    start_test_data_load = time.time()
    X_test, y_test = g.getTensors(X_test_labels, y_test_labels)
    end_test_data_load = time.time()

    callbacks = run_model.createCallbacks()

    start_model_fitting = time.time()
    model = g.createModel(hparams[config.HP_OPTIMIZER], hparams=True)
    log.info('model fitting started')
    #look into ideal batch size
    mod_history = model.fit(X_train, y_train, epochs=10, verbose=True, batch_size=hparams[config.HP_BATCH_SIZE], callbacks = callbacks, validation_split=0.2)
    end_model_fitting = time.time()

    log.info ('model fitting finished successfully')

    timing_measures = [
        end_train_test_split - start_train_test_split,
        end_train_data_load - start_train_data_load,
        end_test_data_load - start_test_data_load,
        end_model_fitting - start_model_fitting,
        end_model_fitting - start_train_test_split
    ]

    log.info('model evaluation started')
    with open('hp_results.txt', 'a') as res_log:
        results = model.evaluate(X_test, y_test, verbose=True)
        res_log.write(' '.join([str(r) for r in results]) + ' \n')
        res_log.write('Timing Benchmarks:\n')
        res_log.write(' '.join([str(r) for r in timing_measures]) + '\n')

    print(results)
    log.info('model evaluation completed')
    #returns loss value (MeanSquaredLogarithmicError) and metric values, currently LogCoshError and coeff_determination
    #and RootMeanSquaredError and MeanSquaredError
    #run_model.createPlot(mod_history)
    return results[0]


def hpRunModel(run_dir, hparams):
  with tf.summary.create_file_writer(run_dir).as_default():
    hp.hparams(hparams)  # record the values used in this trial
    accuracy = hpBuildModel(hparams)
    tf.summary.scalar(config.METRIC_ACCURACY, accuracy, step=1)


def optimizeHyperparameters():
    session_num = 0

#    for num_units in config.HP_NUM_UNITS.domain.values:
#        for dropout_rate in (config.HP_DROPOUT.domain.min_value, config.HP_DROPOUT.domain.max_value):
    for batch_size in config.HP_BATCH_SIZE.domain.values:
        for optimizer in config.HP_OPTIMIZER.domain.values:
            hparams = {
                #config.HP_NUM_UNITS: num_units,
                #config.HP_DROPOUT: dropout_rate,
                config.HP_BATCH_SIZE: batch_size,
                config.HP_OPTIMIZER: optimizer,
            }
            run_name = "run-%d" % session_num
            print('--- Starting trial: %s' % run_name)
            print({h.name: hparams[h] for h in hparams})
            hpRunModel('logs/hparam_tuning/' + run_name, hparams)
            session_num += 1

#optimizeHyperparameters()