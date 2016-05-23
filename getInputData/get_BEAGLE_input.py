#!/usr/bin/env python
import os
import sys
import cProfile
import re

os.chdir("/rsrch1/epi/home/yliu13/AGP_data/")
os.getcwd()

markerDir = "/rsrch1/epi/home/yliu13/AGP_data/MarkerFile_mapping_result/"
dataDir = "/rsrch1/epi/home/yliu13/AGP_data/AGP_original_data_link/"
outputDir = "/rsrch1/epi/home/yliu13/AGP_data/phasedKid_input/"
chrId = sys.argv[1]

# separate 22 chromosomes and organize as BEAGLE input
# from a bunch of "Final_Report_1000_1158.txt"s to each chromosome "BEAGLE_input_chr10.txt"
def getMarkerFile(chrId):
	markerFile = {}
	markerFilename = markerDir+"Marker_chr"+chrId+".txt"
#	id = 0
	for line in open(markerFilename).readlines():
#		id+=1
		snpName,mapInfo,flag,g1,g2 = line.split()
#		markerFile[snpName] = [mapInfo,id]
		markerFile[snpName] = int(mapInfo)
	print "marker file created, size=", len(markerFile)
	return markerFile

def getMarkerList(chrId):
	markerList = []
	markerFilename = markerDir+"Marker_chr"+chrId+".txt"
	for line in open(markerFilename).readlines():
		snpName,mapInfo,flag,g1,g2 = line.split()
		markerList.append(snpName)
	print "marker list created, size=", len(markerList)
	return markerList

def getSampleName(sampleFile):
	famId,sampleId = re.search('_(\d+)_(\d+)\.',sampleFile).groups()
	sampleName = famId+'_'+sampleId
	return sampleName

def getSampleData(filename,markerList):
	sampleData = {}
	sampleName = getSampleName(filename)
	for line in open(filename).readlines():
		if sampleName not in line:
			continue
		item = line.split()
		if item[2] != chrId:
			continue
		snpName = item[1]
		if snpName in markerList:
			sampleData[snpName] = (item[14],item[15])
	print "sample %s data created, size= %d" % (sampleName, len(sampleData))	
	return sampleData,sampleName

def isOrdered(l):
        x = all(l[i] <= l[i+1] for i in xrange(len(l)-1))
        print x

# create a dictionary to save the data in each chromosome and combine them together for BEAGLEinput format
outputDict = {}

# markerData = getMarkerFile(chrId)
markerList = getMarkerList(chrId)

fileList = "kidNameList.txt"
for file in open(fileList).readlines():
	file = file.strip()
	fileName = dataDir + file
	sampleData,sampleName = getSampleData(fileName,markerList)
	outputDict[sampleName] = sampleData

print "ready to write the file"
fw = open(outputDir+"Beagle_input_chr_"+chrId+".txt",'w')
fw.write('I id')
for sampleName in outputDict.keys():
	fw.write(' '+sampleName+' '+sampleName)
fw.write('\n')

for snpName in markerList:
#	print markerData[snpName]
	fw.write('M '+snpName+' ')
	for sampleName in outputDict.keys():
		sampleData = outputDict[sampleName]
		fw.write(' '.join(sampleData[snpName])+' ')
	fw.write('\n')
fw.close()

print "done!"
