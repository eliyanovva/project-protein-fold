import sys
from contextlib import redirect_stdout
sys.path.append('../')
import time
import logging as log

import tensorflow as tf
from tensorboard.plugins.hparams import api as hp

from graph_cnn.model import GraphCNN
import graph_cnn.run_model as run_model
import config


with tf.summary.create_file_writer('logs/hparam_tuning').as_default():
    hp.hparams_config(
        hparams=[
            config.HP_BATCH_SIZE,
            config.HP_DROPOUT,
            config.HP_OPTIMIZER,
            config.HP_LEARNINGRATE,
            config.HP_LOSS,
            config.HP_VALIDATION_SPLIT,
            config.HP_TEST_TRAIN_SPLIT
        ],
        metrics=[hp.Metric(config.METRIC_ACCURACY, display_name='Accuracy')],
    )

class hp_GraphCNN(GraphCNN):
    
    def hp_createModel(self, hparams={
        #config.HP_OPTIMIZER: 'adagrad',
        #config.HP_BATCH_SIZE: 32,
        #config.HP_DROPOUT: 0.15,
        #config.HP_LEARNINGRATE: 0.001,
        config.HP_VALIDATION_SPLIT: 0.15,
        config.HP_TEST_TRAIN_SPLIT: 0.15,
        }):
                
        prot_adj_in = tf.keras.layers.Input(
            shape=(config.PROTEIN_ADJACENCY_MAT_SIZE, config.PROTEIN_ADJACENCY_MAT_SIZE),
            name='Protein-Adjacency-Matrix'
        )
        
        prot_feat_in = tf.keras.layers.Input(
            shape=(config.PROTEIN_ADJACENCY_MAT_SIZE, config.PROTEIN_FEATURES_COUNT),
            name='Protein-Feature-Matrix'
        )

        ligand_adj_in = tf.keras.layers.Input(
            shape=(config.LIGAND_ADJACENCY_MAT_SIZE, config.LIGAND_ADJACENCY_MAT_SIZE),
            name='Ligand-Adjacency-Matrix'
        )

        ligand_feat_in = tf.keras.layers.Input(
            shape=(config.LIGAND_ADJACENCY_MAT_SIZE, config.LIGAND_FEATURES_COUNT),
            name='Ligand-Feature-Matrix'
        )

        dlayer = tf.keras.layers.Dropout(hparams[config.HP_DROPOUT])

        x = tf.keras.layers.Conv1D(filters=1024, kernel_size=3, activation='relu')(prot_adj_in)
        x = tf.keras.layers.MaxPooling1D(pool_size=(2))(x)
        x = tf.keras.layers.Conv1D(filters=512, kernel_size=3, activation='relu')(x)
        x = tf.keras.layers.MaxPooling1D(pool_size=(2))(x)
        x = tf.keras.layers.Conv1D(filters=256, kernel_size=3, activation='relu')(x)
        x = tf.keras.layers.MaxPooling1D(pool_size=(2))(x)
        x = tf.keras.layers.Flatten()(x)
        x = tf.keras.layers.Dense(1024, activation="relu")(x)
        x = dlayer(inputs=x, training=True)
        x = tf.keras.layers.Dense(512, activation="relu")(x)
        x = tf.keras.Model(inputs=prot_adj_in, outputs=x)
        
        y = tf.keras.layers.Flatten()(prot_feat_in)
        y = tf.keras.layers.Dense(512, activation="relu")(y)
        y = dlayer(inputs=y, training=True)
        y = tf.keras.layers.Dense(64, activation="relu")(y)
        y = tf.keras.Model(inputs=prot_feat_in, outputs=y)

        z = tf.keras.layers.Flatten()(ligand_adj_in)
        z = tf.keras.layers.Dense(64, activation="relu")(z)
        z = dlayer(inputs=z, training=True)
        z = tf.keras.layers.Dense(16, activation="relu")(z)
        z = tf.keras.Model(inputs=ligand_adj_in, outputs=z)
        
        z1 = tf.keras.layers.Flatten()(ligand_feat_in)
        z1 = tf.keras.layers.Dense(256, activation="relu")(z1)
        z1 = dlayer(inputs=z1, training=True)
        z1 = tf.keras.layers.Dense(64, activation="relu")(z1)
        z1 = tf.keras.Model(inputs=ligand_feat_in, outputs=z1)

        combined = tf.keras.layers.concatenate([x.output, y.output, z.output, z1.output])
        
        out = tf.keras.layers.Dense(1024, activation="relu")(combined)
        out = tf.keras.layers.Dense(512, activation="relu")(out)
        out = tf.keras.layers.Dense(64, activation="relu")(out)
        out_regression = tf.keras.layers.Dense(1, activation="linear")(out)
        #FIXME: Add the classification layers
        #out_classification = tf.keras.layers.Dense(1, activation='softmax')(out)
        
        model = tf.keras.Model(
            inputs=[x.input, y.input, z.input, z1.input],
            outputs= out_regression#, out_classification]
        )
        
        with open('results.txt', 'a') as res_log:
            with redirect_stdout(res_log):
                model.summary()
            res_log.write('\n')

        model.compile(
            optimizer= hparams[config.HP_OPTIMIZER],
            loss=tf.keras.losses.MeanSquaredError(),
            metrics=[tf.keras.metrics.LogCoshError(),
                tf.keras.metrics.RootMeanSquaredError(),
                ]
        )
        return model

def hpBuildModel(params, hparams={
        config.HP_OPTIMIZER: 'adagrad',
        config.HP_BATCH_SIZE: 32,
        config.HP_DROPOUT: 0.15,
        config.HP_LEARNINGRATE: 0.001,
        config.HP_VALIDATION_SPLIT: 0.15,
        config.HP_TEST_TRAIN_SPLIT: 0.15,
        }):
    g = hp_GraphCNN()
    g.initialize()

    start_train_test_split = time.time()
    X_train_labels, X_test_labels, y_train_labels, y_test_labels = g.trainTestSplit(model_test_size=hparams[config.HP_TEST_TRAIN_SPLIT])
    end_train_test_split = time.time()

    start_train_data_load = time.time()
    X_train, y_train = g.getTensors(X_train_labels, y_train_labels)
    end_train_data_load = time.time()

    start_test_data_load = time.time()
    X_test, y_test = g.getTensors(X_test_labels, y_test_labels)
    end_test_data_load = time.time()

    callbacks = run_model.createCallbacks()

    start_model_fitting = time.time()
    model = g.hp_createModel(hparams)
    log.info('model fitting started')
    mod_history = model.fit(X_train,
        y_train,
        epochs=10,
        verbose=True,
        batch_size=hparams[config.HP_BATCH_SIZE],
        callbacks = callbacks,
        validation_split=hparams[config.HP_VALIDATION_SPLIT])
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
    #explanation = eli5.explain_weights(estimator)
    with open('hp_results.txt', 'a') as res_log:
        results = model.evaluate(X_test, y_test, verbose=True)
        res_log.write('trial for '.join(params) + '\n')
        res_log.write(' '.join([str(r) for r in results]) + ' \n')
        res_log.write('Timing Benchmarks:\n')
        res_log.write(' '.join([str(r) for r in timing_measures]) + '\n')
        #text = eli5.format_as_text(explanation, show_feature_values=True)
        #res_log.write('Format Explanation:\n')
        #res_log.write(' '.join([str(r) for r in text]) + '\n')


    print(results)
    log.info('model evaluation completed')
    #prints loss value (MeanSquaredLogarithmicError) and metric values, currently LogCoshError and coeff_determination
    #and RootMeanSquaredError and MeanSquaredError
    return results[0]
    #returns loss value


def hpRunModel(run_dir, hparams, params):
  with tf.summary.create_file_writer(run_dir).as_default():
    hp.hparams(hparams)  # record the values used in this trial
    accuracy = hpBuildModel(params, hparams)
    tf.summary.scalar(config.METRIC_ACCURACY, accuracy, step=1)


def optimizeHyperparameters():
    session_num = 0

    #for optimizer in config.HP_OPTIMIZER.domain.values:
    #    for batch_size in config.HP_BATCH_SIZE.domain.values:
    #        for dropout in config.HP_DROPOUT.domain.values:
    hparams = {
        config.HP_OPTIMIZER: 'adamax',
            config.HP_LEARNINGRATE: 0.001,
            config.HP_BATCH_SIZE: 256,
            config.HP_DROPOUT: 0.2,
            config.HP_TEST_TRAIN_SPLIT: 0.15,
            config.HP_VALIDATION_SPLIT: 0.15,
        }
    run_name = "run-%d" % session_num
    print('--- Starting trial: %s' % run_name)
    parameters = {h.name: hparams[h] for h in hparams}
    print(parameters)
    hpRunModel('logs/hparam_tuning/' + run_name, hparams, parameters)
    session_num += 1