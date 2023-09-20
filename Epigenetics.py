import math

"This files adds the epigenetic information to the sequence input"


def addEpi(epigentics, ndarray, position):
    for i in range(len(epigentics)):
        if epigentics[i] != 0:
            epi = str(epigentics[i])
            for z in range(len(epi)):
                ndarray[i][z][position] = float(epi[z])
    return ndarray


def addMethylation(epigentics, ndarray, position):
    for i in range(len(epigentics)):
        epi = epigentics[i].split(',')
        floatEpi = [float(e) for e in epi]
        ndarray[i, :, position] = floatEpi
    return ndarray


# For cell types with no methylation information we add the average of T-cell methylation information
def avg_Methylation(epigentics, ndarray, epiLen, position):
    avgf = 0
    nucleotides = len(epigentics[0].split(','))
    for i in range(len(epigentics)):
        epi = epigentics[i].split(',')
        avgf += math.fsum(list(map(float, epi)))
    avgf /= (len(epigentics) * nucleotides)
    for j in range(epiLen):
        METH = [avgf] * nucleotides
        ndarray[j, :, position] = METH
    return ndarray
