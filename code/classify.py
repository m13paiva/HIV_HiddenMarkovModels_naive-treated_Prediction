import numpy as np
import math
import csv

def predicted_probability(pred_score, all_scores):
    # 1) Find max log-score for stable shift
    max_score = max(all_scores)

    # 2) Sum of exponentials of (score - max_score)
    shifted_exps = [np.exp(s - max_score) for s in all_scores]
    denominator = sum(shifted_exps)

    # 3) Numerator is exp(pred_score - max_score)
    numerator = np.exp(pred_score - max_score)

    # 4) Probability
    return numerator / denominator

def normalized_score(pred_clas, all_scores):
    other_clas=max(all_scores)
    epsilon=np.finfo(float).eps

    if pred_clas == other_clas:
        res = 1
    elif pred_clas > 0:
        if other_clas > 0:
            res = other_clas / pred_clas
        else:
            res = epsilon / (pred_clas - other_clas)
    elif pred_clas < 0:
        res = pred_clas / other_clas
    else:
        res = -epsilon / other_clas

    return res ** 4


def classify_seq(seq_dict):
    max_score=-math.inf
    for clas in seq_dict:
        if seq_dict[clas]>=max_score:
            max_score=seq_dict[clas]
            pred_clas=clas
    return pred_clas

def save_output(nested_dict,filename):
    """
        Save a dictionary of dictionaries to a CSV file.

        Parameters:
            nested_dict (dict): The outer dictionary with items as rows, inner dictionary keys as columns.
            filename (str): The filename of the CSV file to save.

        Example Input:
            nested_dict = {
                'row1': {'col1': 1, 'col2': 2},
                'row2': {'col1': 3, 'col2': 4},
            }

        Output CSV:
            col1,col2,row
            1,2,row1
            3,4,row2
        """
    # Ensure all inner dictionaries have the same keys
    all_columns =["seqs"]
    all_columns+=[key for key in next(iter(nested_dict.values())).keys()]

    with open(filename, mode='w', newline='', encoding='utf-8') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=all_columns)

        # Write the header row
        writer.writeheader()

        # Write each row
        for row_key, inner_dict in nested_dict.items():
            row = {**inner_dict, 'seqs': row_key}
            writer.writerow(row)

    print(f"\nOutput data saved successfully to:\n{filename}\n")


def print_output(data):
    for seq in data:
        #print(f"\n»{seq} classified as {data[seq]['pred_class']}, with a predicted probability of {round(data[seq]['pred_prob']*100,2)}%")
        print(f"\n»{seq} classified as {data[seq]['pred_class']}, with a normalized score of {round(data[seq]['norm_score'], 4)}")


def classify(seq_scores):
    seq_classified={}
    for seq in seq_scores:
        seq_classified[seq]={}
        seq_dict=seq_scores[seq]
        pred_clas=classify_seq(seq_dict)
        others={k:seq_dict[k] for k in seq_dict if k!=pred_clas}
        #pred_prob=predicted_probability(seq_dict[pred_clas],list(others.values()))
        norm_score=normalized_score(seq_dict[pred_clas], list(others.values()))
        seq_classified[seq]["pred_class"] = pred_clas
        #seq_classified[seq]["pred_prob"] = pred_prob
        seq_classified[seq]["norm_score"] = norm_score
        seq_classified[seq]={**seq_classified[seq],**seq_dict}
    return seq_classified





