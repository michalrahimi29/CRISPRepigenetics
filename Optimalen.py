import csv
import random
import tensorflow as tf
from CRISPRepi import *

"""This file finds the permutation of flanking sequences that gives the best result"""


# File maker of the 5 cross-validation of every flanking sequences permutation
def optimalLength():
    f = open('optimal Length.csv', 'w', newline='')
    writer = csv.writer(f)
    writer.writerow(['Up', 'Down', 'result'])
    for i in range(14):  # up
        print("******************************" + str(i) + "***********************************")
        for j in range(25):  # down
            acc = findLen(i, j)
            writer.writerow([i, j, acc])


def findLen(up, down):
    data = createTrainSet(up, down, seqs_protospacer, seqs_down, seqs_up, seqs_pam, {})
    k = 5
    fold = int(len(data) / k) + 1
    sumspearman = 0
    for i in range(0, len(data) - 1, fold):
        data_test = data[i:i + fold + 1]
        labels_test = labels[i:i + fold + 1]
        data_train = np.append(data[0:i], data[1 + i + fold:], axis=0)
        weights_train = np.append(weights[0:i], weights[1 + i + fold:], axis=0)
        labels_train = np.append(labels[0:i], labels[1 + i + fold:], axis=0)
        data_train, labels_train = shuffle(data_train, labels_train)
        model1 = leenay_Model(data)
        model1.fit(data_train, labels_train, epochs=10, batch_size=16, verbose=1, sample_weight=weights_train)
        pred_test1 = model1.predict(data_test)
        pred_test = pred_test1
        sumspearman += spearmanr(labels_test, pred_test.reshape(len(pred_test)))[0]
    return sumspearman / k


# Get the permutation of flanking sequences with best result
def getBest():
    a = pd.read_csv("optimal Length.csv")
    vals = a["result"].tolist()
    up = a["Up"].tolist()
    down = a["Down"].tolist()
    maxV = max(vals)
    jindex = vals.index(maxV)
    upLen = up[jindex]
    downLen = down[jindex]
    print(str(upLen) + "," + str(downLen) + "," + str(maxV))


if __name__ == '__main__':
    a = pd.read_csv("Final_leenay_dataset.csv")
    seqs_protospacer = a["protospacer"]
    seqs_pam = a["PAM"].tolist()
    seqs_up = a["upstream"].tolist()
    seqs_down = a["downstream"].tolist()
    no_var = a["no_variant"].to_numpy()
    reads = a["total_reads"].to_numpy()
    Total = np.sum(reads)
    W = np.divide(reads, Total)
    epsilon = 0.2
    weights = W + epsilon
    l = np.divide(no_var, reads)
    labels = np.subtract(1.0, l)
    #optimalLength()
    getBest()
