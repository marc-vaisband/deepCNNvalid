import numpy as np
import sys
import os
import keras
import tensorflow
import random
from keras.utils import to_categorical
sys.path.append(os.path.abspath(".."))
from utils.scoring_utils import score_from_cat
from training.make_keras_model import make_batchnorm_model

"""
I/O
"""
print("Loading data.")
contexted_data_folder = os.path.abspath("../data/in_facility/contexted")
contexted_tensor_path = os.path.join(contexted_data_folder, "scrambled_x1k2_tensor.npy")
contexted_labels_path = os.path.join(contexted_data_folder, "x1k2_labels.npy")
with open(contexted_tensor_path, "rb") as f:
    train_data = np.load(f)
with open(contexted_labels_path, "rb") as f:
    train_labels_num = np.load(f)

train_data = keras.utils.normalize(train_data)
train_labels = to_categorical(train_labels_num)

validation_data_path = os.path.abspath("../data/kotani/scrambled_kotani_x1k2_tensor.npy")
validation_label_path = os.path.abspath("../data/kotani/kotani_labels.npy")
with open(validation_data_path, "rb") as f:
    valid_data = np.load(f)
with open(validation_label_path, "rb") as f:
    valid_labels_num = np.load(f)

valid_labels = to_categorical(valid_labels_num)
valid_data = keras.utils.normalize(valid_data)


"""
First: Train model on all in-facility data. 
"""
print("Training model:")
fit_args = {"x": train_data,
            "y": train_labels,
            "batch_size": 256,
            "epochs": 50,
            "verbose": 1,
            }

np.random.seed(0)
tensorflow.random.set_seed(0)
random.seed(0)


model = make_batchnorm_model(input_shape=train_data.shape[1:], metrics=["accuracy"])
history = model.fit(**fit_args)

model_out_folder = os.path.abspath("../misc/keras_model_for_kotani")
os.makedirs(model_out_folder, exist_ok=True)
model.save(model_out_folder)

"""
Next: Evaluate predictions on external data.
"""
y_pred_cat = model.predict(valid_data)
print("Final validation metrics:")
acc, precision, recall, f1 = score_from_cat(y_true_cat=valid_labels, y_pred_cat=y_pred_cat, verbose=True)

report_folder = os.path.abspath("../reports/kotani_validation")
os.makedirs(report_folder, exist_ok=True)
with open(os.path.join(report_folder, "validation_metrics.csv"), "w") as f:
    f.write(";".join(["accuracy", "precision", "recall", "f1"]) + "\n")
    f.write(f"{acc};{precision};{recall};{f1}\n")

