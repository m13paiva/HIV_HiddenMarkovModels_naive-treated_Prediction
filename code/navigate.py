import os
import pandas as pd
import csv
from pathlib import Path
from align import Align

def get_abs_path(path1,path2=None,path3=None,path4=None):
    # Get the absolute path of the current script's directory
    current_dir = Path(__file__).resolve().parent

    # Navigate to the main folder (parent of the code folder)
    main_folder = current_dir.parent

    if path2:
        if path3:
            if path4:
                return main_folder / path1 / path2 / path3 /path4
            else:
                return main_folder / path1 / path2 / path3
        else:
            return main_folder / path1 / path2

    else:
        return main_folder / path1

def save_dict_to_txt_file(data, file_path):
    """
    Saves a dictionary to a .txt file.

    Args:
        data (dict): The dictionary to save.
        file_path (str): The path of the .txt file to save to.
    """
    try:
        with open(file_path, 'w') as file:
            for key, value in data.items():
                file.write(f"{key}: {value}\n")
    except Exception as e:
        print(f"An error occurred: {e}")

def read_txt_file_to_dict(file_path):
    """
    Reads a .txt file where each line is in the format "key: value" or "key: [value1, value2, ...]"
    and returns a dictionary with those key-value pairs. List values are parsed into Python lists.

    Args:
        file_path (str): The path of the .txt file to read from.

    Returns:
        dict: A dictionary containing the key-value pairs read from the file.
    """
    import ast  # For safely evaluating list-like strings
    data_dict = {}
    try:
        with open(file_path, 'r') as file:
            for line in file:
                # Strip whitespace from the line
                line = line.strip()

                # Skip any empty lines
                if not line:
                    continue

                # Split only on the first colon to allow colons in value
                parts = line.split(':', 1)

                # If the line doesn't contain a colon, skip it
                if len(parts) < 2:
                    continue

                key, value = parts[0].strip(), parts[1].strip()

                # Check if the value looks like a list
                if value.startswith('[') and value.endswith(']'):
                    try:
                        # Safely evaluate the list-like string
                        parsed_value = ast.literal_eval(value)
                        if isinstance(parsed_value, list):
                            data_dict[key] = parsed_value
                        else:
                            data_dict[key] = value  # Store as-is if not a list
                    except Exception:
                        data_dict[key] = value  # Store as-is if parsing fails
                else:
                    # Store as a regular string
                    data_dict[key] = value
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")

    return data_dict


def read_fasta(file_path, dis_drm=False, codons=None):
    """
    Reads a FASTA file and returns a list of tuples (header, sequence).
    Each header is a string (minus the leading '>'), and the sequence
    is a single concatenated string of all lines following the header.

    Processing steps on each sequence:
      1. If dis_drm is True, remove characters at positions listed in `codons`
         (1-based) from the raw sequence (which may still contain dashes).
      2. Remove all '-' characters.
      3. Prepend '-' to the resulting sequence.

    Args:
        file_path (str): The path to the FASTA file.
        dis_drm (bool, optional): Whether to exclude specific positions. Defaults to False.
        codons (list, optional): A list of integer positions (1-based) to remove if dis_drm is True.

    Returns:
        list of tuples: A list where each element is (header, sequence).
    """
    if codons is None:
        codons = []

    sequences_list = []
    current_label = None
    current_seq_lines = []

    try:
        with open(file_path, "r") as f:
            for line in f:
                line = line.strip()
                # Skip empty lines
                if not line:
                    continue

                # Check if line is a header (starts with '>')
                if line.startswith(">"):
                    # If there's a label already being tracked, store the previous sequence
                    if current_label is not None:
                        # Combine all lines for the current sequence
                        joined_seq = "".join(current_seq_lines)

                        # (1) If dis_drm is True, remove positions listed in codons (1-based)
                        if dis_drm and codons:
                            to_remove = set(pos - 1 for pos in codons)  # convert to 0-based
                            joined_seq = "".join(
                                char for i, char in enumerate(joined_seq)
                                if i not in to_remove
                            )

                        # (2) Remove all '-' characters
                        joined_seq = joined_seq.replace("-", "")

                        # (3) Prepend '-'
                        joined_seq = "-" + joined_seq

                        # Add to our list of sequences
                        sequences_list.append((current_label, joined_seq))

                    # Start tracking a new header/sequence
                    current_label = line[1:].strip()  # remove the '>' and extra whitespace
                    current_seq_lines = []
                else:
                    # It's part of a sequence
                    current_seq_lines.append(line)

            # After the loop, handle the last sequence if it exists
            if current_label is not None:
                joined_seq = "".join(current_seq_lines)

                # (1) If dis_drm is True, remove positions listed in codons (1-based)
                if dis_drm and codons:
                    to_remove = set(pos - 1 for pos in codons)  # convert to 0-based
                    joined_seq = "".join(
                        char for i, char in enumerate(joined_seq)
                        if i not in to_remove
                    )

                # (2) Remove all '-' characters
                joined_seq = joined_seq.replace("-", "")

                # (3) Prepend '-'
                joined_seq = "-" + joined_seq

                sequences_list.append((current_label, joined_seq))

    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")

    return sequences_list

def save_choices(d):
    file_path = get_abs_path("model_data", "choices.txt")
    save_dict_to_txt_file(d, file_path)

def delete_choices():
    file_path = get_abs_path("model_data", "choices.txt")
    os.remove(file_path)


def access_choices():
    file_path = get_abs_path("model_data", "choices.txt")
    return read_txt_file_to_dict(file_path)

def read_msa(file_path):
    """The myAlign will have the list of sequences from the msa_file and
      the family name is given by the name of the file. The created object is
      then added to the list of families (self.families)."""
    with open(file_path, 'r') as file:
        file_name = os.path.basename(file_path)  # Extracts the file name cross-platform
        subtype = file_name.split('.')[0]  # Splits by '.' and gets the first part
        seq_list = []
        for line in file:
            clean_line = line.strip()  # Remove whitespace, including newlines
            if clean_line:  # Only add non-empty lines
                seq_list.append(clean_line)
        msa = Align(seq_list, subtype)
        return msa

def read_test(folder_path, seq_count):
    for file_path in os.listdir(folder_path):
        with open(os.path.join(folder_path,file_path), 'r') as file:
            file_name = os.path.basename(file_path)  # Extracts the file name cross-platform
            subtype = file_name.split('.')[0]  # Splits by '.' and gets the first part
            seq_list = []
            for line in file:
                clean_line = line.strip()  # Remove whitespace, including newlines
                if clean_line:  # Only add non-empty lines
                    seq_list.append((f"seq_{seq_count}",subtype,"-"+clean_line))
                    seq_count+=1

        return seq_list, seq_count

def filter_max_score(d):
    max_key = max(d, key=d.get)
    max_value = d[max_key]

    return (max_key, max_value)

def write_output(profiles, scores, file_name):
    # Verificar se o arquivo já existe
    if os.path.exists(file_name):
        print(f"The file '{file_name}' already exists and will be replaced.")

    rows = []
    first_row = ["subtype", "pred_subtype", "correct_pred"]
    for subtype in profiles:
        first_row.append(subtype)
    rows.append(first_row)

    for seq in scores:
        row = []
        subtype = seq[1]
        pred_subtype = filter_max_score(scores[seq])[0]
        row.append(subtype)
        row.append(pred_subtype)
        correct_pred = 1 if subtype == pred_subtype else 0
        row.append(correct_pred)
        for score in scores[seq].values():
            row.append(score)
        rows.append(row)

    # Escrever os dados no arquivo CSV
    with open(file_name, mode='w', newline='', encoding='utf-8') as file:
        writer = csv.writer(file)
        writer.writerows(rows)
        #print(f"File '{file_name}' created with success.")

def read_output(file):
    df=pd.read_csv(file)
    l = list(df["correct_pred"])
    accuracy = l.count(1) / len(l)
    subtypes=[subtype for subtype in df.columns[3:]]
    performance={"precision":{},"recall":{},"f1-score":{}}
    for subtype in subtypes:
        positive_df=df[df["pred_subtype"]==subtype]
        negative_df=df[df["pred_subtype"]!=subtype]
        true_pos=len(positive_df[positive_df["subtype"]==subtype])
        false_pos = len(positive_df[positive_df["subtype"] != subtype])
        false_neg=len(negative_df[negative_df["subtype"]==subtype])
        denom=true_pos+false_pos
        if denom!=0:
            precision=true_pos/denom
        else:
            precision=0
        denom=true_pos + false_neg
        if denom!=0:
            recall=true_pos / denom
        else:
            recall=0
        performance["precision"][subtype] = precision
        performance["recall"][subtype] = recall
        denom=precision+recall
        if denom!=0:
            performance["f1-score"][subtype] = 2*(precision*recall/denom)
        else:
            performance["f1-score"][subtype] = 0

    for param in performance:
        performance[param]["mean"]=sum(performance[param].values())/len(performance[param])
    print(f"Accuracy: {accuracy*100}%")
    print(pd.DataFrame(performance))








