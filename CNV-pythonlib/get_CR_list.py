#!/usr/bin/env python
import os
import sys
import numpy
import myLib

fileListName = sys.argv[1]
missingLabel = sys.argv[2]
outputFilename = sys.argv[3]

fileList = [x.strip() for x in open(fileListName).readlines()]

fw = open(outputFilename,'w')
print('\t'.join(['sampleName','GT_call_rate']),file=fw)
for filename in fileList:
        missingCR = myLib.calMissingCallRate(filename, missingLabel)
        print('\t'.join([os.path.basename(filename),str(missingCR)]),file=fw)

fw.close()
print("done")
