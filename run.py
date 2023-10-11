from CRISPRepi import *
import tensorflow as tf
import random
import sys

"This file responsible for initialization and running"


def initialize():
    randomNumSeed = 123
    np.random.seed(randomNumSeed)
    tf.random.set_seed(0)
    random.seed(10)
    epigenetics = sys.argv[1:]
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
    epigeneticDic = {}  # dictionary of epigenetic information, every marker has a specific key
    for i in range(len(epigenetics)):
        b = pd.read_csv("epigenetics_" + epigenetics[i] + ".csv")
        if epigenetics[i] == 'chromatin_accessibility':
            name = 'epi1'
        elif epigenetics[i] == 'CTCF_binding':
            name = 'epi2'
        elif epigenetics[i] == 'H3K4ME3':
            name = 'epi3'
        elif epigenetics[i] == 'methylation':
            name = 'epi4'
        epigeneticDic[name] = b["epigenetics"].to_numpy()
    #lennaysRun(seqs_protospacer, seqs_pam, seqs_up, seqs_down, weights, labels, epigeneticDic)
    lennayPredicionOnHumanCells("hek293", seqs_protospacer, seqs_down, seqs_up, seqs_pam, labels, weights, epigeneticDic)


def save_trained_model(seqs_protospacer, seqs_pam, seqs_up, seqs_down, weights, labels, epigeneticDic):
    data = createTrainSet(seqs_protospacer, seqs_down, seqs_up, seqs_pam, epigeneticDic)
    modeli = leenay_Model(data)
    modeli.fit(data, labels, epochs=25, batch_size=16, verbose=1, sample_weight=weights)
    modeli.save("CRISPRepi_model.keras")
    #reconstructed_model = keras.models.load_model("CRISPRepi_model.keras")  # use this line for uploading the trained model


if __name__ == '__main__':
    initialize()
