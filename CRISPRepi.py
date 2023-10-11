import pandas as pd
import numpy as np
from keras.models import Sequential
from keras.layers import *
from scipy.stats.stats import pearsonr
from scipy.stats.stats import spearmanr
from sklearn.utils import shuffle
from Epigenetics import *
import keras

SIZE = 8  # size of matrix input
UP = 0  # optimal permutation
DOWN = 9


def oneHot(string):
    trantab = str.maketrans('ACGT', '0123')
    string = str(string)
    data = [int(x) for x in list(string.translate(trantab))]
    ret = np.eye(SIZE)[data]
    return ret


def leenay_Model(data):
    model1 = Sequential()
    model1.add(Flatten(input_shape=(data.shape[1], data.shape[2])))
    model1.add(Dense(64, activation='sigmoid'))
    model1.add(Dense(8, activation='sigmoid'))
    model1.add(Dense(1, activation='sigmoid'))
    model1.compile(optimizer='adam', loss='binary_crossentropy')
    return model1


# creates the data for training and evaluating the model
def createTrainSet(seqs_protospacer, seqs_down, seqs_up, seqs_pam, epigeneticDic):
    data_protospacer = np.array(list(map(oneHot, seqs_protospacer)))
    downSeqs = []
    upSeqs = []
    for seq in seqs_down:
        downSeqs.append(seq[:DOWN])
    for seq in seqs_up:
        upSeqs.append(seq[:UP])
    data_pam = np.array(list(map(oneHot, seqs_pam)))
    data_down = np.array(list(map(oneHot, np.array(downSeqs))))
    data_up = np.array(list(map(oneHot, np.array(upSeqs))))
    data = np.append(data_up, data_protospacer, axis=1)
    data = np.append(data, data_pam, axis=1)
    data = np.append(data, data_down, axis=1)
    index = 4
    # adding epigenetic information to sequence input
    for key in epigeneticDic.keys():
        if key == 'epi4':  # specific key for methylation
            data = addMethylation(epigeneticDic[key], data, index)
        else:
            data = addEpi(epigeneticDic[key], data, index)
        index += 1
    return data


# creates the test set of other cell types for evaluation of the model
def humanTestCell(file, cell_type, epigeneticDic):
    a = pd.read_csv(file)
    temp = a["sequence"].tolist()
    seqs = []
    for elem in temp:
        seqs.append(elem[30:62])
    data = np.array(list(map(oneHot, seqs)))
    index = 4
    # adding epigenetic information to sequence input
    for key in epigeneticDic.keys():
        if key == 'epi1':
            b = pd.read_csv(cell_type + "_chromatin_accessibility.csv")
            Newepi = b["epigenetics"].to_numpy()
            data = addEpi(Newepi, data, index)
            index += 1
        elif key == 'epi2':
            c = pd.read_csv(cell_type + "_CTCF_binding.csv")
            Newepi2 = c["epigenetics"].to_numpy()
            data = addEpi(Newepi2, data, index)
            index += 1
        elif key == 'epi3':
            d = pd.read_csv(cell_type + "_H3K4ME3.csv")
            Newepi3 = d["epigenetics"].to_numpy()
            data = addEpi(Newepi3, data, index)
            index += 1
        elif key == 'epi4':
            data = avg_Methylation(epigeneticDic['epi4'], data, len(Newepi), index)
            index += 1
    return data


# 5-Cross validation for training and evaluation
def lennaysRun(protospacer, pam, up, down, weights, labels, epigenetics):
    data = createTrainSet(protospacer, down, up, pam, epigenetics)
    k = 5
    fold = int(len(data) / k) + 1
    sumpearson = 0
    sumspearman = 0
    for i in range(1, len(data) - 1, fold):
        data_test = data[i:i + fold]
        labels_test = labels[i:i + fold]
        data_train = np.append(data[1:i], data[1 + i + fold:], axis=0)
        weights_train = np.append(weights[1:i], weights[1 + i + fold:], axis=0)
        labels_train = np.append(labels[1:i], labels[1 + i + fold:], axis=0)
        data_train, labels_train = shuffle(data_train, labels_train)
        model1 = leenay_Model(data)
        model1.fit(data_train, labels_train, epochs=10, batch_size=16, verbose=1, sample_weight=weights_train)
        pred_test1 = model1.predict(data_test)
        pred_test = pred_test1
        sumpearson = sumpearson + pearsonr(labels_test, pred_test.reshape(len(pred_test)))[0]
        sumspearman = sumspearman + spearmanr(labels_test, pred_test.reshape(len(pred_test)))[0]
        print(pearsonr(labels_test, pred_test.reshape(len(pred_test))))
    print(sumspearman / k)


# used for evaluation of the model on different cell types. works on k562+hek293+hct116+H1
def lennayPredicionOnHumanCells(cell_type, epigenetics):
    data_test = humanTestCell("hek293+k562+hct116_efficicency.csv", cell_type, epigenetics)
    a = pd.read_csv("hek293+k562+hct116_efficicency.csv")
    key = "efficiency_" + cell_type
    labels_test = a[key].to_numpy()
    """
    # in case some guides don't have efficiency value for specific cell_type
    to_delete = np.argwhere(np.isnan(labels_test)).flatten()
    labels_test = np.delete(labels_test, to_delete)
    data_test = np.delete(data_test, to_delete, 0)
    data_train, labels_train = shuffle(data, labels)
    """
    model1 = keras.models.load_model("CRISPRepi_model.keras")
    pred_test1 = model1.predict(data_test)
    pred_test = pred_test1
    print(spearmanr(labels_test, pred_test.reshape(len(pred_test))))
