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
print('\t'.join(['sampleName','BAF_std','LRR_std']),file=fw)
for filename in fileList:
        BAF_sd = myLib.calBAF_stat(filename,missingLabel)[2]
        LRR_sd = myLib.calLRR_stat(filename,missingLabel)[2]
        print('\t'.join([os.path.basename(filename),str(BAF_sd),str(LRR_sd)]),file=fw)

fw.close()
print("done")
