import csv
import statistics
from _csv import reader

OptLEN = 32

"This file parse bigwig files and creates an epigenetic methylation file"


def parseBigWig():
    dictionary = {}
    for i in range(1, 4):
        with open(f"methyl per base individual{i}.txt") as f:
            lines = f.readlines()
            for j in range(len(lines)):
                lineList = lines[j].split('\t')
                chr = lineList[0][3:]
                start = int(lineList[1])
                end = int(lineList[2])
                score = float(lineList[3])
                if chr in dictionary.keys():
                    if (start, end) in dictionary[chr].keys():
                        if score != -1.0:
                            num1, num2 = dictionary[chr][(start, end)]
                            dictionary[chr][(start, end)] = (num1 + score, num2 + 1)
                    else:
                        if score != -1.0:
                            dictionary[chr][(start, end)] = (score, 1)
                        else:
                            dictionary[chr][(start, end)] = (0, 0)
                else:
                    if score != -1.0:
                        dictionary[chr] = {(start, end): (score, 1)}
                    else:
                        dictionary[chr] = {(start, end): (0, 0)}
    return dictionary


def methylationFile(dict):
    vals = []
    for chr in dict.keys():
        for interval in dict[chr].keys():
            if dict[chr][interval][1] != 0:
                avg = dict[chr][interval][0] / dict[chr][interval][1]
                vals.append(avg)
            else:
                avg = -1.0
            dict[chr][interval] = avg
    median = statistics.median(vals)
    for chr in dict.keys():
        for interval in dict[chr].keys():
            if dict[chr][interval] == -1.0:
                dict[chr][interval] = median
    return dict


def createEpi(methylation):
    epi = {}
    with open("Final_leenay_dataset.csv", 'r') as read_obj:
        csv_reader = reader(read_obj)
        next(csv_reader)
        line = 1
        for row in csv_reader:
            chromNum = row[0]
            startG = int(row[1]) + 13
            endG = int(row[2]) - 16
            temp = ["0"] * OptLEN
            for tup in methylation[chromNum].keys():
                startM = tup[0]
                endM = tup[1]
                # the start of the guide
                if startM <= startG < endM <= endG:
                    for i in range(endM - startG + 1):
                        temp[i] = str(methylation[chromNum][tup])
                # the middle of guide
                elif startG <= startM < endM <= endG:
                    for i in range(endM - startM + 1):
                        temp[startM - startG + i] = str(methylation[chromNum][tup])
                # the end of guide
                elif startG <= startM < endG <= endM:
                    for i in range(endG - startM + 1):
                        temp[startM - startG + i] = str(methylation[chromNum][tup])
                    break
                elif startM > endG:
                    break
            epi[line] = ','.join(temp)
            line += 1
    with open('epigenetics_methylation.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['name', 'epigenetics'])
        for k in epi.keys():
            writer.writerow([k, epi[k]])


if __name__ == '__main__':
    dictionary = parseBigWig()
    updated = methylationFile(dictionary)
    print("yay")
    createEpi(updated)