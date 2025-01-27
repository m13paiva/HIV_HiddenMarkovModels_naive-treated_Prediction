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

def test():
    from HMM_prediction.read_raw_data import process
    data=process()
    train, test = generate_train_test(data,0.3)
    print("»»»»»»»»»»»»»»»»»»»»»TRAIN»»»»»»»»»»»»»»»»»»»»")
    print(train)
    print("»»»»»»»»»»»»»»»»»»»»»»TEST»»»»»»»»»»»»»»»»»»»»")
    print(test)
#test()