import numpy as np

def score_from_num(y_true_num, y_pred_num, verbose=False):
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
    assert y_true_cat.shape[0] == y_pred_cat.shape[0]
    assert y_true_cat.shape[1] == y_pred_cat.shape[1] == 2
 
    y_true = np.argmax(y_true_cat, axis=1)
    y_pred = np.argmax(y_pred_cat, axis=1)

    return score_from_num(y_true_num=y_true, y_pred_num=y_pred, verbose=verbose)

