from sklearn.model_selection import StratifiedKFold
import numpy as np
import keras
import tensorflow
import sys
import os
import time
import datetime
import inspect
from make_keras_model import make_batchnorm_model
sys.path.append(os.path.abspath(".."))
from utils.scoring_utils import score_from_cat



contexted_data_folder = os.path.abspath("../data/in_facility/contexted")
contexted_tensor_path = os.path.join(contexted_data_folder, "scrambled_x1k2_tensor.npy")
contexted_labels_path = os.path.join(contexted_data_folder, "x1k2_labels.npy")

nocontext_data_folder = os.path.abspath("../data/in_facility/nocontext")
nocontext_tensor_path = os.path.join(nocontext_data_folder, "scrambled_nocontext_tensor.npy")
nocontext_labels_path = os.path.join(nocontext_data_folder, "nocontext_labels.npy")


def run_5fold_cv(data_path, label_path, get_model_callable, model_name=""):
    """
    Function that performs stratified k-fold crossvalidation

    :param data_path: Path to binary dump of numpy tensor
    :param label_path: Path to binary dump of labels
    :param get_model_callable: Callable that takes the arguments 'input_shape' and 'metrics' and returns a compiled
    keras model ready for fitting.
    :param model_name: Name of folder in which reports and models will be stored
    :return:
    """

    # I/O
    execution_datetime = datetime.datetime.fromtimestamp(time.time())

    concrete_out_folder = \
        os.path.join(os.path.abspath("../stored_models/cv/"),
                     f"{model_name}_"
                     f"{execution_datetime.year}_{execution_datetime.month}_{execution_datetime.day}_"
                     f"{execution_datetime.hour}_{execution_datetime.minute}_{execution_datetime.second}/")
    model_save_folder = os.path.join(concrete_out_folder, "models")
    info_save_folder = os.path.join(concrete_out_folder, "info")
    os.makedirs(concrete_out_folder, exist_ok=True)
    os.makedirs(model_save_folder, exist_ok=True)
    os.makedirs(info_save_folder, exist_ok=True)

    print("Loading data.")
    with open(data_path, "rb") as f:
        in_data = np.load(f)

    with open(label_path, "rb") as f:
        labels = np.array(np.load(f))

    cat_labels = keras.utils.to_categorical(labels)



    skf_seed = 42
    np.random.seed(skf_seed)
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=skf_seed)

    valaccs = []
    train_seeds = []
    precisions = []
    recalls = []
    f1s = []

    i = 0
    for train_idc, test_idc in skf.split(in_data, labels):
        train_seed = i + 42
        np.random.seed(train_seed)
        tensorflow.random.set_seed(train_seed)
        train_seeds.append(train_seed)

        # Train-test-split
        train_data, test_data = keras.utils.normalize(in_data[train_idc, ...]), keras.utils.normalize(
            in_data[test_idc, ...])

        train_labels, test_labels = cat_labels[train_idc], cat_labels[test_idc]

        train_metrics = ["accuracy"]

        model = get_model_callable(input_shape=train_data.shape[1:], metrics=train_metrics)

        fit_args = {"x": train_data,
                    "y": train_labels,
                    "batch_size": 256,
                    "epochs": 50,
                    "verbose": 0,
                    }

        history = model.fit(**fit_args)

        model.save(os.path.join(model_save_folder, str(i)))

        model_prediction_test = model.predict(test_data)

        valacc, precision, recall, f1 = score_from_cat(y_true_cat=test_labels, y_pred_cat=model_prediction_test)

        valaccs.append(valacc)
        precisions.append(precision)
        recalls.append(recall)
        f1s.append(f1)

        print(f"Fold{i + 1}. Accuracy: {valacc}, precision: {precision}, recall: {recall}, f1: {f1}")
        i += 1

    print(valaccs)
    print(f"Average validation accuracy: {np.average(valaccs)}")
    print(f"Std of valaccs: {np.std(valaccs)}")
    print("Saving settings and results")

    with open(info_save_folder + "paths.txt", "w") as f:
        f.write(f"Data path: {data_path}\n")
        f.write(f"Labels path: {label_path}\n")

    with open(info_save_folder + "code.txt", "w") as f:
        f.write(inspect.getsource(get_model_callable))

    with open(info_save_folder + "seeds.txt", "w") as f:
        f.write(f"Skf seed: {skf_seed}\n")
        train_seed_string = ";".join(list(map(str, train_seeds)))
        f.write(f"Train seeds: {train_seed_string}")

    with open(info_save_folder + "args.txt", "w") as f:
        for key in ["batch_size", "epochs"]:
            f.write(f"{key}: {fit_args[key]}\n")

    with open(info_save_folder + "metrics.txt", "w") as f:
        f.write("metric;scores;avg\n")
        f.write("accuracy;" + ",".join([str(val) for val in valaccs]) + f";{np.mean(valaccs)}\n")
        f.write("precision;" + ",".join([str(val) for val in precisions]) + f";{np.mean(precisions)}\n")
        f.write("recall;" + ",".join([str(val) for val in recalls]) + f";{np.mean(recalls)}\n")
        f.write("f1;" + ",".join([str(val) for val in f1s]) + f";{np.mean(f1s)}\n")


run_5fold_cv(data_path=nocontext_tensor_path, label_path=nocontext_labels_path,
             get_model_callable=make_batchnorm_model, model_name="nocontext")
run_5fold_cv(data_path=contexted_tensor_path, label_path=contexted_labels_path,
             get_model_callable=make_batchnorm_model, model_name="contexted")

