import numpy as np

def score_from_num(y_true_num, y_pred_num, verbose=False):
    """
    Utility function that takes binary (!) true and predicted labels and converts them into accuracy, precision,
    recall and F1 score, to evaluate prediction quality.

    :param y_true_num: np.array of true numerical class labels
    :param y_pred_num: np.array of predicted numerical class labels
    :param verbose: bool to signify whether the scores should be printed to console
    :return: Accuracy, Precision, Recall and F1 score calculated
    """


    assert y_true_num.shape[0] == y_pred_num.shape[0]
    tp_mask = np.logical_and(y_pred_num == 1, y_true_num == 1)
    fp_mask = np.logical_and(y_pred_num == 1, y_true_num == 0)
    tn_mask = np.logical_and(y_pred_num == 0, y_true_num == 0)
    fn_mask = np.logical_and(y_pred_num == 0, y_true_num == 1)

    tp = np.sum(tp_mask)
    fp = np.sum(fp_mask)
    fn = np.sum(fn_mask)
    tn = np.sum(tn_mask)

    acc = (tp + tn) / len(y_true_num)
    precision = tp / (tp + fp)
    recall = tp / (tp + fn)
    f1 = 2 * (precision * recall) / (precision + recall)

    if verbose:
        print(f"TP: {tp}, FP: {fp}, TN: {tn}, FN: {fn}")
        print(f"Acc {acc:.2f}, precision {precision:.2f}, recall {recall:.2f}, f1 {f1:.2f}")

    return acc, precision, recall, f1


def score_from_cat(y_true_cat, y_pred_cat, verbose=False):
    """
    Utility function that takes binary (!) true and predicted labels in categorical form and converts them into
    accuracy, precision, recall and F1 score, to evaluate prediction quality.

    :param y_true_cat: np.array of true categorical class labels
    :param y_pred_cat: np.array of predicted categorical class labels
    :param verbose: bool to signify whether the scores should be printed to console
    :return: Accuracy, Precision, Recall and F1 score calculated
    """

    assert y_true_cat.shape[0] == y_pred_cat.shape[0]
    assert y_true_cat.shape[1] == y_pred_cat.shape[1] == 2
 
    y_true = np.argmax(y_true_cat, axis=1)
    y_pred = np.argmax(y_pred_cat, axis=1)

    return score_from_num(y_true_num=y_true, y_pred_num=y_pred, verbose=verbose)

