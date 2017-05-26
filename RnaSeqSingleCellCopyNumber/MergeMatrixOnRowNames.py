# Author: Stephen Piccolo

import os, sys, glob
import utilities

inFilePath1 = sys.argv[1]
inFilePath2 = sys.argv[2]
outFilePath = sys.argv[3]

data1 = utilities.readMatrixFromFile(inFilePath1)
data2 = utilities.readMatrixFromFile(inFilePath2)

header1 = data1.pop(0)
header2 = data2.pop(0)
if len(header2) == len(data2[0]):
    header2.pop(0)

headerCombined = header1 + header2

data1Dict = {}
for row in data1:
    data1Dict[row[0]] = data1Dict.setdefault(row[0], []) + [row]

data2Dict = {}
for row in data2:
    data2Dict[row[0]] = data2Dict.setdefault(row[0], []) + [row[1:]]

outFile = open(outFilePath, 'w')
outFile.write("\t".join(headerCombined) + "\n")
for key in sorted(data1Dict.keys()):
    if data1Dict.has_key(key) and data2Dict.has_key(key):
        for row1 in data1Dict[key]:
            for row2 in data2Dict[key]:
                outFile.write("\t".join(row1 + row2) + "\n")
outFile.close()
