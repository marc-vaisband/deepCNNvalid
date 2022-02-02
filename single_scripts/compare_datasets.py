import numpy as np
import sys
import os
import keras
from keras.models import load_model as tensorflow_load
import pickle
from collections import Counter
from keras.utils import to_categorical, normalize
sys.path.append(os.path.abspath("/limcr-ngs/Marc/git/PhD/limcr_mutcall/"))
from mutcallClean.utils import scoring_utils
from mutcallClean.training.model_building import nol2_models
from sklearn.model_selection import train_test_split

sys.path.append(os.path.abspath("/limcr-ngs/Marc/git/PhD/limcr_mutcall/mutcallClean/data_generation"))
import pairLocus

# our_data_path = "/limcr-ngs/Marc/git/PhD/limcr_mutcall/mutcallClean/data_generation/data/x1k2_depth_noindels/x1k2_noindels_depth_tensor.npy"
# our_label_path = "/limcr-ngs/Marc/git/PhD/limcr_mutcall/mutcallClean/data_generation/data/x1k2_depth_noindels/x1k2_noindels_depth_labels.npy"
# with open(our_data_path, "rb") as f:
#     train_data = np.load(f)
# with open(our_label_path, "rb") as f:
#     train_labels_num = np.load(f)

# train_labels = to_categorical(train_labels_num)

comparison_data_path = "/limcr-ngs/Marc/git/PhD/limcr_mutcall/mutcallClean/data_generation/data/testset1/nonsyn_depth_tensor.npy"
comparison_label_path = "/limcr-ngs/Marc/git/PhD/limcr_mutcall/mutcallClean/data_generation/data/testset1/nonsyn_labels.npy"
comparison_ploci_path = "/limcr-ngs/Marc/git/PhD/limcr_mutcall/mutcallClean/data_generation/data/testset1/nonsyn_ploci.pkl"
with open(comparison_data_path, "rb") as f:
    valid_data = np.load(f)
with open(comparison_label_path, "rb") as f:
    valid_labels_num = np.load(f)
with open(comparison_ploci_path, "rb") as f:
    valid_ploci = pickle.load(f)

valid_labels = to_categorical(valid_labels_num)


"""
# Uncomment this to kick the comparison bams
depth_cutoff = train_data.shape[-1] // 2
assert depth_cutoff == 14

train_data = train_data[..., :depth_cutoff]
valid_data = valid_data[..., :depth_cutoff]
"""




# train_data = keras.utils.normalize(train_data)
valid_data = keras.utils.normalize(valid_data)

# get_model_callable = nol2_models.flatten_batchnorm_model
# fit_args = {"x": train_data,
#             "y": train_labels,
#             "batch_size": 256,
#             "epochs": 50,
#             # "validation_data": (valid_data, valid_labels),
#             "verbose": 1,
#             }
# train_metrics = ["accuracy"]




"""
First: Train model on all of our data, then see. 
"""

# model = get_model_callable(input_shape=train_data.shape[1:], metrics=train_metrics)
# history = model.fit(**fit_args)
#
model_out_folder = "/limcr-ngs/Marc/git/PhD/limcr_mutcall/mutcallClean/misc/model_on_our_data"
# os.makedirs(model_out_folder, exist_ok=True)
# model.save(model_out_folder)

model = tensorflow_load(model_out_folder)

y_pred_cat = model.predict(valid_data)
print("Final valid metrics:")
_ = scoring_utils.score_from_cat(y_true_cat=valid_labels, y_pred_cat = y_pred_cat, verbose=True)

y_pred_num = np.argmax(y_pred_cat, axis=1)


# with open("/limcr-ngs/Marc/git/PhD/limcr_mutcall/mutcallClean/misc/testset1_all_mutations_according_to_NN.csv", "w") as f:
#     f.write("Mouse;GL_ID;CL_ID;chr;start;end;ref;alt;predicted_mut_prob\n")
#     for current_pl, mut_prob in zip(valid_ploci, y_pred_cat[:, 1]):
#         if mut_prob > 0.5:
#             current_mouse = current_pl.kwdict["mouse_ID"]
#             f.write(f"{current_mouse};{current_pl.GL_ID};{current_pl.CL_ID};{current_pl.chromosome};{current_pl.start};{current_pl.stop};{current_pl.ref};{current_pl.alt};{mut_prob}\n")
# #

# with open("/limcr-ngs/Marc/git/PhD/limcr_mutcall/mutcallClean/misc/testset1_all_fns.csv", "w") as f:
#     f.write("Mouse;GL_ID;CL_ID;chr;start;end;ref;alt;predicted_mut_prob\n")
#     for current_pl, label, mut_prob in zip(valid_ploci, valid_labels_num, y_pred_cat[:, 1]):
#         if mut_prob < 0.5 and label == 1:
#             current_mouse = current_pl.kwdict["mouse_ID"]
#             f.write(f"{current_mouse};{current_pl.GL_ID};{current_pl.CL_ID};{current_pl.chromosome};{current_pl.start};{current_pl.stop};{current_pl.ref};{current_pl.alt};{mut_prob}\n")
#

fp_ploci = []

with open("/limcr-ngs/Marc/git/PhD/limcr_mutcall/mutcallClean/misc/testset1_nonsyn_fps.csv", "w") as f:
    f.write("Mouse;GL_ID;CL_ID;chr;start;end;ref;alt;gene;predicted_mut_prob\n")
    for current_pl, mut_prob, label in zip(valid_ploci, y_pred_cat[:, 1], valid_labels_num):
        if mut_prob > 0.5 and label == 0:
            fp_ploci.append(current_pl)
            current_mouse = current_pl.kwdict["mouse_ID"]
            f.write(f"{current_mouse};{current_pl.GL_ID};{current_pl.CL_ID};{current_pl.chromosome};{current_pl.start};{current_pl.stop};{current_pl.ref};{current_pl.alt};{current_pl.gene};{mut_prob}\n")

print("Exonic func counter:")
print(Counter([pl.kwdict["exonic_func"] for pl in fp_ploci]))

print("Gene detail counter:")
print(Counter([pl.gene_detail for pl in fp_ploci]))
