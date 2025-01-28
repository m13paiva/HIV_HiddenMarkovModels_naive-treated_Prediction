import os
import pandas as pd
from sklearn.model_selection import train_test_split

def split(df,test_fraction,random_state=None):
    train, test = train_test_split(df, test_size=test_fraction, random_state=random_state)
    return train, test

def generate_train_test(data,test_fraction,random_state=None):
    train={}
    test={}
    for clas in data:
        train[clas]={}
        test[clas]={}
        for subtype in data[clas]:
            train[clas][subtype], test[clas][subtype] = split(data[clas][subtype],test_fraction,random_state)
    return train, test

def read_and_clean_csv(file):
    df = pd.read_csv(file, sep=';', index_col=None, header=0)
    df.columns = list(df.columns[:3]) + df.columns[3:].astype(int).tolist()
    df = df.dropna()
    df = df.iloc[:, 1:]

    return df


def remove_ambiguity(df):
    positions = df[df == "?"].stack().index.tolist()
    for pos in positions:
        df.at[pos[0], pos[1]] = df.iloc[0, pos[1] + 1]  # pos[1]+1 because iloc works with indexes (starting from 0)
    return df.drop(0, axis=0)


def remove_major_mutations(df,codons=None):
    #print(codons)
    return df.drop(codons, axis=1)


def split_treatment(df):
    df_naive = df[df["Treatment"] == "Naive"]
    df_treated = df[df["Treatment"] == "Treated"]
    df_naive = df_naive.drop("Treatment", axis=1)
    df_treated = df_treated.drop("Treatment", axis=1)
    return df_naive, df_treated


def split_by_subtypes(df):
    d = {}
    subtypes = df["Subtype"].unique()
    for subtype in subtypes:
        sub_df = df[df["Subtype"] == subtype].drop("Subtype", axis=1)
        d[subtype] = sub_df
    return d

def equilibrate_dataset(data,rs=None):
    for subtype in data:
        min_len=min(len(df) for df in data[subtype].values())
        for clas in data[subtype]:
            df=data[subtype][clas]
            if len(df)>min_len:
                data[subtype][clas]=df.sample(n=min_len,random_state=rs)

def process(codons,file="global_AA_table_wts.csv", remove_major_muts=True, split_subtypes=True, eq=True, rs=None):
    df = read_and_clean_csv(file)
    df = remove_ambiguity(df)
    if remove_major_muts:
        df = remove_major_mutations(df,codons)
    if split_subtypes:
        d=split_by_subtypes(df)
    else:
        d={"unknown":df.drop("Subtype", axis=1)}
    data={}
    for subtype in d:
        data[subtype]={}
        df_naive, df_treated = split_treatment(d[subtype])
        data[subtype]["naive"]=df_naive
        data[subtype]["treated"] = df_treated
    if eq:
        equilibrate_dataset(data,rs)
    return data




def check_all_ambiguity_cols(df):
    cols = []
    for i, col in df.itercol():
        for char in col:
            if char != "?":
                break
        cols.append(i)
    return cols

def downsize_df(df,n_samples, rs=None):
    len_df=len(df)
    if n_samples>len_df:
        n_samples=len_df
    return df.sample(n=n_samples,random_state=rs)

def downsize_data(data,n_samples, random_state=None):
    for clas in data:
        for subtype in data[clas]:
            data[clas][subtype]=downsize_df(data[clas][subtype],
                                            n_samples,
                                            random_state)
    return data

def df_to_str_list(df,rem_gaps=False):
    l=[]
    for _,row in df.iterrows():
        str_l=[]
        if rem_gaps:
            str_l=[char for char in row if char!="-"]
        else:
            str_l=[char for char in row]

        delimiter = "" # Define a delimiter
        join_str = delimiter.join(str_l)
        l.append(join_str)
    return l


def df_to_txt(df, file_path, rem_gaps=False):
    """
    Write the contents of a DataFrame to a text file. If the file already exists, it will be replaced.

    Args:
        df (DataFrame): The DataFrame to be converted to text.
        file_path (str): The path where the text file will be saved.
        rem_gaps (bool): If True, removes gaps from the text.
    """
    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(file_path), exist_ok=True)

    # Convert DataFrame to string list
    str_list = df_to_str_list(df, rem_gaps)

    # Write to the file
    with open(file_path, 'w', encoding='utf-8') as f:
        for item in str_list:
            f.write(item + '\n')


def save_data(data,folder,subtype,rem_gaps=False,n_samples=None,random_state=None):
    # Delete all files in the folder
    if os.path.exists(folder):
        for file in os.listdir(folder):
            file_path = os.path.join(folder, file)
            if os.path.isfile(file_path):
                os.remove(file_path)
    if n_samples:
        #print(n_samples)
        if subtype!="unknown":
            for clas in data[subtype]:
                df=downsize_df(data[subtype][clas],n_samples,random_state)
                #print(f"LENNN {len(df)}")
                df_to_txt(df,f"{folder}/{subtype}_{clas}.txt",rem_gaps)
        else:
            for subtype_ in data:
                for clas in data[subtype_]:
                    df = downsize_df(data[subtype_][clas], n_samples, random_state)
                    df_to_txt(df,f"{folder}/{subtype_}_{clas}.txt",rem_gaps)
    else:
        if subtype!="unknown":
            for clas in data[subtype]:
                df_to_txt(data[subtype][clas],f"{folder}/{subtype}_{clas}.txt",rem_gaps)
        else:
            for subtype_ in data:
                for clas in data[subtype_]:
                    df_to_txt(data[subtype_][clas],f"{folder}/{subtype_}_{clas}.txt",rem_gaps)

def prepare_hmm_data(rawdata_path,proccessdata_path,codons, subtype,
                     remove_major_mutations=True, split_subtypes=True, eq=True, rs=None):
    data = process(codons,rawdata_path, remove_major_mutations, split_subtypes, eq, rs)
    save_data(data, proccessdata_path, subtype, rem_gaps=False, n_samples=None)

def prepare_hmm_data_test(file_path,train_path,test_path,codons, subtype,
                          test_fraction, n_test_samples=None, n_train_samples=None,
                          remove_major_mutations=True,split_subtypes=True,eq=True,rs=None):
    #print(f"test {n_test_samples} train {n_train_samples}")
    data = process(codons,file_path, remove_major_mutations, split_subtypes, eq, rs)
    train, test = generate_train_test(data, test_fraction,rs)
    save_data(train, train_path, subtype, rem_gaps=False, n_samples=n_train_samples,random_state=rs)
    save_data(test, test_path, subtype, rem_gaps=True, n_samples=n_test_samples,random_state=rs)

