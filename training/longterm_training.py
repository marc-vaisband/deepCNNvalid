import numpy as np
import keras
import tensorflow
import sys
import os
import time
import datetime
import inspect
from sklearn.model_selection import train_test_split
from make_keras_model import make_batchnorm_model
sys.path.append(os.path.abspath(".."))
from utils.scoring_utils import score_from_cat


contexted_data_folder = os.path.abspath("../data/in_facility/contexted")
contexted_tensor_path = os.path.join(contexted_data_folder, "scrambled_x1k2_tensor.npy")
contexted_labels_path = os.path.join(contexted_data_folder, "x1k2_labels.npy")

nocontext_data_folder = os.path.abspath("../data/in_facility/nocontext")
nocontext_tensor_path = os.path.join(nocontext_data_folder, "scrambled_nocontext_tensor.npy")
nocontext_labels_path = os.path.join(nocontext_data_folder, "nocontext_labels.npy")



def run_longterm_evaluation(data_path, label_path, get_model_callable, model_name="", n_splits=100, tf_seed=0, skf_seed=0):
    """
    Function that performs evaluation by averaging performance metrics from a large number of random
    stratified train-test-splits.

    :param data_path: Path to binary dump of numpy tensor
    :param label_path: Path to binary dump of labels
    :param get_model_callable: Callable that takes the arguments 'input_shape' and 'metrics' and returns a compiled
    keras model ready for fitting.
    :param model_name: Name of folder in which reports and models will be stored
    :param n_splits: Number of random stratified train-test-splits to perform
    :param tf_seed: Seed for tensorflow
    :param skf_seed: Seed for K-fold
    :return:
    """

    execution_datetime = datetime.datetime.fromtimestamp(time.time())

    concrete_out_folder = \
        os.path.join(os.path.abspath("../reports/100splits/"),
                     f"{model_name}_"
                     f"{execution_datetime.year}_{execution_datetime.month}_{execution_datetime.day}_"
                     f"{execution_datetime.hour}_{execution_datetime.minute}_{execution_datetime.second}/")

    info_save_folder = concrete_out_folder
    os.makedirs(concrete_out_folder, exist_ok=True)
    os.makedirs(info_save_folder, exist_ok=True)

    print("Loading data.")
    with open(data_path, "rb") as f:
        in_data = np.load(f)

    with open(label_path, "rb") as f:
        labels = np.array(np.load(f))

    cat_labels = keras.utils.to_categorical(labels)



    valaccs = []
    precisions = []
    recalls = []
    f1s = []

    split_seeds = []

    np_tf_seed = 0
    np.random.seed(np_tf_seed)
    tensorflow.random.set_seed(np_tf_seed)

    for i in range(n_splits):
        split_seed = i
        split_seeds.append(split_seed)

        train_idc, test_idc = train_test_split(np.arange(len(in_data)), train_size=0.8, shuffle=True,
                                               random_state=i, stratify=labels)

        train_data = keras.utils.normalize(in_data[train_idc, ...])
        test_data = keras.utils.normalize(in_data[test_idc, ...])

        train_labels, test_labels = cat_labels[train_idc], cat_labels[test_idc]

        test_classes = np.argmax(test_labels, axis=1)
        assert len(test_classes) == len(test_labels)

        train_metrics = ["accuracy"]

        model = get_model_callable(input_shape=train_data.shape[1:], metrics=train_metrics)

        fit_args = {"x": train_data,
                    "y": train_labels,
                    "batch_size": 256,
                    "epochs": 50,
                    "verbose": 0,
                    }
        np.random.seed(split_seed)
        tensorflow.random.set_seed(split_seed)
        history = model.fit(**fit_args)

        model_prediction_test = model.predict(test_data)

        valacc, precision, recall, f1 = score_from_cat(y_true_cat=test_labels, y_pred_cat=model_prediction_test)

        valaccs.append(valacc)
        precisions.append(precision)
        recalls.append(recall)
        f1s.append(f1)

        i += 1
        print(f"Done with split number {i}. Current average valacc: {np.mean(valaccs)}")

    print(valaccs)
    print(f"Average validation accuracy: {np.average(valaccs)}")
    print(f"Std of valaccs: {np.std(valaccs)}")
    print("Saving settings and results")

    with open(os.path.join(info_save_folder,  "code.txt"), "w") as f:
        f.write(inspect.getsource(get_model_callable))

    with open(os.path.join(info_save_folder,  "seeds.txt"), "w") as f:
        f.write(f"Numpy/Tensorflow seed: {np_tf_seed}\n")
        split_seed_string = ";".join(list(map(str, split_seeds)))
        f.write(f"Train seeds: {split_seed_string}")

    with open(os.path.join(info_save_folder,  "args.txt"), "w") as f:
        for key in ["batch_size", "epochs"]:
            f.write(f"{key}: {fit_args[key]}\n")

    with open(os.path.join(info_save_folder,  "metrics.txt"), "w") as f:
        f.write("metric;avg;std;scores\n")
        f.write("accuracy;" f"{np.mean(valaccs)};" + f"{np.std(valaccs)};" + ",".join([str(val) for val in valaccs]) + "\n")
        f.write("precision;" + f"{np.mean(precisions)};" + f"{np.std(precisions)};" + ",".join([str(val) for val in precisions]) + "\n")
        f.write("recall;" + f"{np.mean(recalls)};" + f"{np.std(recalls)};" + ",".join([str(val) for val in recalls]) + "\n")
        f.write("f1;" + f"{np.mean(f1s)};" + f"{np.std(f1s)};" + ",".join([str(val) for val in f1s]) + "\n")

    print(valaccs)
    print(f"Average validation accuracy: {np.average(valaccs)}")
    print(f"Std of valaccs: {np.std(valaccs)}")


run_longterm_evaluation(data_path=contexted_tensor_path, label_path=contexted_labels_path,
                        get_model_callable=make_batchnorm_model, model_name="contexted")
run_longterm_evaluation(data_path=nocontext_tensor_path, label_path=nocontext_labels_path,
                        get_model_callable=make_batchnorm_model, model_name="nocontext")


