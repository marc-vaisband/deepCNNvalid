import keras
import keras.layers as layers
import tensorflow as tf

def make_batchnorm_model(input_shape, metrics=None):
    if metrics is None:
        metrics = ["accuracy"]

    kernelsize_1 = (3, 3)

    kernelsize = kernelsize_1
    model_CNN = keras.models.Sequential(name="flatten_batchnorm")

    model_CNN.add(layers.Conv2D(32, kernel_size=kernelsize, activation='relu', input_shape=input_shape,
                                data_format="channels_last"))
    model_CNN.add(layers.MaxPooling2D(pool_size=kernelsize, padding="same"))

    model_CNN.add(layers.Conv2D(16, kernel_size=kernelsize, activation='relu'))
    model_CNN.add(layers.MaxPooling2D(pool_size=kernelsize, padding="same"))

    model_CNN.add(layers.Conv2D(32, kernel_size=kernelsize, activation='relu'))
    model_CNN.add(layers.MaxPooling2D(pool_size=kernelsize, padding="same"))

    model_CNN.add(layers.Conv2D(32, kernel_size=(1, 1), activation='relu'))
    model_CNN.add(layers.BatchNormalization(momentum=0.8))
    model_CNN.add(layers.Flatten())

    model_CNN.add(layers.Dense(10, activation='relu'))
    model_CNN.add(layers.Dropout(0.2))

    model_CNN.add(layers.Dense(2, activation='softmax'))

    model_CNN.compile(loss=tf.keras.losses.BinaryCrossentropy(label_smoothing=0.1), optimizer="adam", metrics=metrics)
    return model_CNN
