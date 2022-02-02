import numpy as np
from sklearn.model_selection import StratifiedKFold
import sys
import os
import time
import random
import keras
import tensorflow
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

x1k3_noindels_data_path = "/limcr-ngs/Marc/git/PhD/limcr_mutcall/mutcallClean/data_generation/data/x1k3_depth_noindels/x1k3_noindels_depth_tensor.npy"
x1k3_noindels_labels_path = "/limcr-ngs/Marc/git/PhD/limcr_mutcall/mutcallClean/data_generation/data/x1k3_depth_noindels/x1k3_noindels_depth_labels.npy"
noindels_ploci_path = "/limcr-ngs/Marc/git/PhD/limcr_mutcall/mutcallClean/data_generation/data/ploci.pkl"

with open(x1k3_noindels_data_path, "rb") as f:
    x1k3_data = np.load(f)
with open(x1k3_noindels_labels_path, "rb") as f:
    num_labels = np.load(f)
with open(noindels_ploci_path, "rb") as f:
    noindels_ploci = pickle.load(f)

cat_labels = to_categorical(num_labels)

flatten_batchnorm_callable = nol2_models.flatten_batchnorm_model

skf_seed = 42
np.random.seed(skf_seed)
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=skf_seed)

valid_predictions = np.zeros(shape=cat_labels.shape)
valaccs = []

i = 0
for train_idc, test_idc in skf.split(x1k3_data, num_labels):
    model_out_folder = f"/limcr-ngs/Marc/git/PhD/limcr_mutcall/mutcallClean/misc/5foldcv_models/{i}_withnorm"
    os.makedirs(model_out_folder, exist_ok=True)

    train_seed = i+42
    random.seed(train_seed)
    np.random.seed(train_seed)
    tensorflow.random.set_seed(train_seed)

    train_data, test_data = x1k3_data[train_idc, ...], x1k3_data[test_idc, ...]
    train_data, test_data = keras.utils.normalize(train_data), keras.utils.normalize(test_data)

    train_labels, test_labels = cat_labels[train_idc], cat_labels[test_idc]

    test_classes = np.argmax(test_labels, axis=1)
    assert len(test_classes) == len(test_labels)
    assert all(test_classes == num_labels[test_idc])



    train_metrics = ["accuracy"]

    model = flatten_batchnorm_callable(input_shape=train_data.shape[1:], metrics=train_metrics)

    fit_args = {"x": train_data,
                "y": train_labels,
                "batch_size": 256,
                "epochs": 50,
                "verbose": 1,
                }
    history = model.fit(**fit_args)

    model.save(model_out_folder)

    # model = tensorflow_load(model_out_folder)

    model_prediction_test = model.predict(test_data)
    predicted_test_classes = np.argmax(model_prediction_test, axis=1)

    valid_predictions[test_idc] = model_prediction_test
    valacc, _, _, _ = scoring_utils.score_from_cat(y_true_cat = test_labels, y_pred_cat = model_prediction_test)
    print(f"Fold {i+1} validation accuracy: {valacc}")
    valaccs.append(valacc)

    i += 1
    



print(valaccs)
print(f"Average: {np.average(valaccs)}")


overall_valacc, _, _, _ = scoring_utils.score_from_cat(y_true_cat = cat_labels, y_pred_cat = valid_predictions)
print(f"For cross-checking: Average accuracy with all validation predictions: {overall_valacc}")


print("Saving predictions.")
with open("/limcr-ngs/Marc/git/PhD/limcr_mutcall/mutcallClean/misc/5foldcv_models/all_valid_predictions_withnorm.npy", "wb") as f:
    np.save(f, valid_predictions)
print("Generating list of mistakes")

wrong_ploci = []

with open("/limcr-ngs/Marc/git/PhD/limcr_mutcall/mutcallClean/misc/5foldcv_models/all_misclassified_withnorm.csv", "w") as f:
    f.write("GL_ID;CL_ID;chr;start;end;ref;alt;ref-alt;funcrefgene;predicted_mut_prob\n")
    for current_pl, mut_prob, label in zip(noindels_ploci, valid_predictions[:, 1], num_labels):
        if mut_prob > 0.5 and label == 0 or mut_prob <= 0.5 and label == 1:
            wrong_ploci.append(current_pl)
            f.write(f"{current_pl.GL_ID};{current_pl.CL_ID};{current_pl.chromosome};{current_pl.start};{current_pl.stop};{current_pl.ref};{current_pl.alt};{current_pl.ref}-{current_pl.alt};{current_pl.funcrefgene};{mut_prob}\n")



print("Initialising idx dicts")
base_strings = ["A", "C", "T", "G"]
pl_alt_idc = {base_string: [] for base_string in base_strings}
pl_ref_idc = {base_string: [] for base_string in base_strings}

funcrefgene_entries = list(set([pl.funcrefgene for pl in noindels_ploci]))
pl_funcrefgene_idc = {frg: [] for frg in funcrefgene_entries}

print("Filling indices")

for k, pl in enumerate(noindels_ploci):
    pl_alt_idc[pl.alt].append(k)
    pl_ref_idc[pl.ref].append(k)
    pl_funcrefgene_idc[pl.funcrefgene].append(k)

for base_string in base_strings:
    pl_alt_idc[base_string] = np.array(pl_alt_idc[base_string])
    pl_ref_idc[base_string] = np.array(pl_ref_idc[base_string])

for frg in pl_funcrefgene_idc.keys():
    pl_funcrefgene_idc[frg] = np.array(pl_funcrefgene_idc[frg])

print("Building accuracy dicts..")
print("Alt, ref")
# Accuracies based on alt
alt_accs = {}
ref_accs = {}
for base_string in base_strings:
    alt_mask = pl_alt_idc[base_string]
    alt_accs[base_string] = scoring_utils.score_from_cat(y_true_cat = cat_labels[alt_mask], y_pred_cat = valid_predictions[alt_mask])[0]
    
    ref_mask = pl_ref_idc[base_string]
    ref_accs[base_string] = scoring_utils.score_from_cat(y_true_cat = cat_labels[ref_mask], y_pred_cat = valid_predictions[ref_mask])[0]



print("Refalt")
ref_alt_idc = {}
ref_alt_accs = {}
for ref in base_strings:
    for alt in base_strings:
        if ref == alt:
            ref_alt_accs[(ref, alt)] = np.nan
            continue

        idc = np.intersect1d(pl_ref_idc[ref], pl_alt_idc[alt])
        ref_alt_idc[(ref, alt)] = idc
        ref_alt_accs[(ref, alt)] = scoring_utils.score_from_cat(y_true_cat = cat_labels[idc], y_pred_cat = valid_predictions[idc])[0]

# Accuracies based on funcrefgene
print("Funcrefgene")
funcrefgene_accs = {frg: scoring_utils.score_from_cat(y_true_cat = cat_labels[pl_funcrefgene_idc[frg]], y_pred_cat = valid_predictions[pl_funcrefgene_idc[frg]])[0] for frg in funcrefgene_entries}

print("Start dumping")
with open("/limcr-ngs/Marc/git/PhD/limcr_mutcall/mutcallClean/misc/category_idx_objects3.pkl", "wb") as f:
    pickle.dump((pl_ref_idc, pl_alt_idc, pl_funcrefgene_idc, ref_alt_idc), f)
with open("/limcr-ngs/Marc/git/PhD/limcr_mutcall/mutcallClean/misc/category_acc_objects3.pkl", "wb") as f:
    pickle.dump((ref_accs, alt_accs, funcrefgene_accs, ref_alt_accs), f)


print("Done dumping")
#

