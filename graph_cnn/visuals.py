import matplotlib.pyplot as plt


def plotLossCurve(mod_history):
    #plot the loss curve: test vs training
    fig, ax1 = plt.subplots(1, figsize=(15, 5))

    ax1.plot(mod_history.history["loss"])
    ax1.plot(mod_history.history["val_loss"])
    ax1.legend(["train", "test"], loc="upper right")
    ax1.set_xlabel("Epochs")
    ax1.set_ylabel("Loss")

    plt.savefig('visuals/loss_graph.png')
    plt.close()

def scatterPlotLosses(mod_history):
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