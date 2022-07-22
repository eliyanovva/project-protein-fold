import time
from graph_cnn.model import GraphCNN
import config
import logging as log
import matplotlib.pyplot as plt
import tensorflow as tf
import visuals


def createPlot(mod_history):
    #plot the loss curve: test vs training
    fig, ax1 = plt.subplots(1, figsize=(15, 5))

    ax1.plot(mod_history.history["loss"])
    ax1.plot(mod_history.history["val_loss"])
    ax1.legend(["train", "test"], loc="upper right")
    ax1.set_xlabel("Epochs")
    ax1.set_ylabel("Loss")

    plt.savefig('visuals/loss_graph.png')
    plt.close()


def createScatterPlot(mod_history):
    # Plot the results
    #print(mod_history.history.keys())
    #acc = mod_history.history['accuracy']
    loss = mod_history.history['loss']
    val_loss = mod_history.history['val_loss']
    epochs = range(len(loss))

    plt.scatter(epochs[1:], loss[1:], c='red', label='Training MSE')
    plt.scatter(epochs[1:], val_loss[1:], c='blue', label='Validation MSE')
    plt.title('Training accuracy')
    plt.legend(loc=0)

    plt.savefig('visuals/training_accuracy.png')
    plt.close()


def createCallbacks():
    callbacks = [
        tf.keras.callbacks.EarlyStopping(
            monitor = 'val_loss',
            patience = 0,
            restore_best_weights = True,
            mode = 'min',
        )
    ]

    return callbacks


def trainGNN(gnn, X_train, y_train, tr_batch_size=32, tr_optimizer='adagrad'):
    model = gnn.createModel(hp_optimizer=tr_optimizer)
    callbacks = createCallbacks()

    log.info('model fitting started')
    mod_history = model.fit(
        X_train, y_train, epochs=10, verbose=True,
        batch_size=tr_batch_size, callbacks = callbacks, validation_split=0.2)
    
    log.info ('model fitting finished successfully')

    return model, mod_history


def testGNN(model, X_test, y_test):
    log.info('model evaluation started')
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


def runModel(batch_size=-1, test_frac=0.3):
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
    model, training_history = trainGNN(gnn, X_train, y_train)
    end_model_fitting = time.time()

    start_model_fitting = time.time()
    results = testGNN(model, X_test, y_test)
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
