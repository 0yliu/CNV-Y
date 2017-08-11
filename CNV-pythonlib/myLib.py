import os
import sys
import numpy as np
import numpy
import pandas
import pandas as pd
import re
import random
from scipy.stats import norm
from pandas import Series,DataFrame
from operator import itemgetter
from itertools import groupby
from datetime import datetime
import dupCalc

## List all the files in a directory
def listFiles(directory):
	# list all the files in DIRECTORY
	files=os.listdir(directory)
	#filter non result files here if necessary
	files = filter(lambda x:x[0]!='.', files)
	return files

def getNumOfMarkers(startSNP,endSNP,chrId,markerFile):
	# for the standard input files, there is no number of markers info.
	# When do simulations, we want to use the same number of markers to observed the distribution
	# get the number of markers info based on the startSNP and endSNP
	# markerFile = 'source/Cran_clean_markers_simpleList.txt'
	#SNP.Name	Chr	MapInfo
	#rs3094315	1	742429
	#rs3131972	1	742584
	df = pandas.read_table(markerFile,index_col=0)
	df['rowNum'] = range(len(df))
	if df.loc[startSNP,'Chr'] != int(chrId):
		print("ERROR: the startSNP is not located in the indicated chromosome.")
		sys.exit(1)
	if df.loc[startSNP,'Chr'] != df.loc[endSNP,'Chr']:
		print("ERROR: the startSNP and the endSNP are not in the same chromosome.")
		sys.exit(1)
	x = df.loc[endSNP,'rowNum'] - df.loc[startSNP,'rowNum']+1
	return(x)

def getOrderedMarkers(markerFile):
	df = pd.read_table(markerFile,index_col=0)
	df['rowNum'] = range(len(df))
	return(df)

## get a dictionary about pedigree structure (more importantly, the trio structure)
def getTrioDict(trioList,dataSource):
	## the dict key is the offspring name;
	## the list corresponding to each key is father,mother,tissue_type
	trioDict = {}
	if dataSource == 'craniofacial':
		for l in open(trioList).readlines():
			if 'filename-offspring' in l:
				continue
			items = l.split()
			s1 = re.split('[\_|\.]',items[1]); offspringName = s1[5]
			s2 = re.split('[\_|\.]',items[4]); fatherName = s2[5]
			s3 = re.split('[\_|\.]',items[7]); motherName = s3[5]
			trioDict[offspringName] = fatherName,motherName,items[2]
	if dataSource == 'autism':
		for l in open(trioList).readlines():
			if 'filename-offspring' in l:
				continue
			items = l.split()
			s1 = re.split('[\_|\.]',items[1]); offspringName = s1[2]+'_'+s1[3]
			s2 = re.split('[\_|\.]',items[4]); fatherName = s2[2]+'_'+s2[3]
			s3 = re.split('[\_|\.]',items[7]); motherName = s3[2]+'_'+s3[3]
			trioDict[offspringName] = fatherName,motherName,items[2]
	return(trioDict)

## extract PennCNV results and ordered within chromosomes
def getBadRegions(pennCNVDir,kidName):
	rmMkrs = {}
	for chrId in range(1,23):
		rmMkrs['chr'+str(chrId)] = {}
	pennCNVFilename = pennCNVDir + kidName + '.adjusted.rawcnv'
	for cnvLine in open(pennCNVFilename).readlines():
		items = cnvLine.strip().split()
		chrNum = re.split(':',items[0])[0]
		tmpList = rmMkrs[chrNum]
		eventSampleName = re.split('\-|\.',items[4])[1]
		rmStartSNP = re.split('=',items[5])[1]
		rmEndSNP = re.split('=',items[6])[1]
		startLoc = re.split('\:|\-',items[0])[1]
		tmpList[int(startLoc)] = [rmStartSNP,rmEndSNP]
		rmMkrs[chrNum] = tmpList
	return rmMkrs

## extract PennCNV results and ordered within chromosomes
def getParentalBadRegions(pennCNVFilename,parentsNames_pair):
	rmMkrs = {}
	for chrId in range(1,23):
		rmMkrs['chr'+str(chrId)] = {}
	for cnvLine in open(pennCNVFilename).readlines():
		items = cnvLine.strip().split()
		chrNum = re.split(':',items[0])[0]
		tmpList = rmMkrs[chrNum]
		eventSampleName = re.split('\-|\.',items[4])[1]
		if eventSampleName in parentsNames_pair:
			rmStartSNP = re.split('=',items[5])[1]
			rmEndSNP = re.split('=',items[6])[1]
			startLoc = re.split('\:|\-',items[0])[1]
			tmpList[int(startLoc)] = [rmStartSNP,rmEndSNP]
		rmMkrs[chrNum] = tmpList
	return rmMkrs

## Extract the chromosome length and then I can select the length of the region and the start SNP position
def getNormalChrRegions(markersDF,chrId,rmMkrs):
	CNVDict = DataFrame(rmMkrs['chr'+chrId],index=['startSNP','endSNP']) # colnames are sorted keys
	subDF = markersDF[markersDF.Chr==int(chrId)]
	## generate normal markers list to sample
	# If the CNVDict is empty, then extract the entire chromosome as available region. 
	if CNVDict.empty:	
		includeMkrsList = [list(subDF.index)]
		#includeMkrsList = [list(subDF['rowNum'])]
	else:
		rmMkrDF = DataFrame(columns=['startLoc','endLoc','startSNP','endSNP'])
		workingLine = 0
		for badStartLoc in CNVDict.columns.values.tolist():
			rmStartSNP,rmEndSNP = CNVDict[badStartLoc]
			rmStartLoc = subDF[subDF['SNP.Name']==rmStartSNP].index[0]
			rmEndLoc = subDF[subDF['SNP.Name']==rmEndSNP].index[0]
			#rmStartLoc = subDF.loc[rmStartSNP,'rowNum']
			#rmEndLoc = subDF.loc[rmEndSNP,'rowNum']
			if workingLine>0 and rmMkrDF.loc[workingLine-1,'endLoc']>=rmStartLoc:
				rmMkrDF.loc[workingLine-1,'endLoc'] = rmEndLoc
				rmMkrDF.loc[workingLine-1,'endSNP'] = rmEndSNP
			else:
				rmMkrDF.loc[workingLine] = rmStartLoc,rmEndLoc,rmStartSNP,rmEndSNP
				workingLine +=1
		rmMkrsList = []
		for idx in rmMkrDF.index:
			rmMkrsList.extend(list(range(int(rmMkrDF.loc[idx,'startLoc']),int(rmMkrDF.loc[idx,'endLoc'])+1)))
		includeMkrsList = []
		#for k,g in groupby(enumerate(set(subDF['rowNum'])-set(rmMkrsList)), lambda x:x[0]-x[1]):
		for k,g in groupby(enumerate(set(subDF.index)-set(rmMkrsList)), lambda x:x[0]-x[1]):
			includeMkrsList.append(list(map(itemgetter(1),g)))
	return(includeMkrsList,subDF)

## get the output information for generating inputData of our dup model
def getNormalRegionOutput(includeMkrsList,subDF,length):
	lenList = [len(x) for x in includeMkrsList]
	passList = list(filter(lambda y: y[1]>length, enumerate(lenList)))
	if len(passList)==0:
		return(None)
	else:
		selectMkrsListIndex = random.choice(passList)[0]
		selectMkrsDF = subDF.loc[includeMkrsList[selectMkrsListIndex]]
		startSNPidx = random.choice(selectMkrsDF.index[:-length])
		endSNPidx = startSNPidx + length - 1
		dupInfo = [selectMkrsDF.loc[startSNPidx,'SNP.Name'],selectMkrsDF.loc[startSNPidx,'MapInfo'],selectMkrsDF.loc[endSNPidx,'SNP.Name'],selectMkrsDF.loc[endSNPidx,'MapInfo']]
	return(dupInfo)

def simulate_gt(dupDF,sd1,sd2):
	## sd1 is the sd for the homogeneous markers 'AAA' and 'BBB'
	## sd2 is the sd for the heterozygous markers 'AAB' and 'ABB'
	kidGT = []
	for snpName in dupDF.index:
		GT = dupDF.ix[snpName,'gt']
		if GT=='AAA':
			x = abs(np.random.normal(0,sd1))
		elif GT=='AAB':
			x = abs(np.random.normal(1/3,sd2))
		elif GT=='ABB':
			x = abs(np.random.normal(2/3,sd2))
		elif GT=='BBB':
			x = abs(np.random.normal(1,sd1))
			if x>1:
				x=2-x
		kidGT.append(x)
	return(kidGT)

def estimate_params(pennCNVDir,gtDataDir,sampleName,cnvName, perc, maxNumsOfMkrs):
	## go into sampleName and get the normal region (at most maxNumsOfMkrs)
	## then extract the homozygous AA and BB markers baf information, as well as AB markers baf info.
	## calculate the distribution parameters
	rmMkrs = getBadRegions(pennCNVDir, cnvName)
	sampleDF = pd.read_table(gtDataDir+sampleName+'.txt')
	includeMkrsList = []
	for chrId in range(1,23):
		CNVDict = DataFrame(rmMkrs['chr'+str(chrId)],index=['startSNP','endSNP'])	
		subDF = sampleDF[sampleDF.Chr==int(chrId)]
		if CNVDict.empty:
			includeMkrsList.extend(list(subDF.index))
		else:
			workingLine = 0
			rmMkrDF = DataFrame(columns=['startLoc','endLoc','startSNP','endSNP'])
			for badStartLoc in CNVDict.columns.values.tolist():
				rmStartSNP,rmEndSNP = CNVDict[badStartLoc]
				rmStartLoc = subDF[subDF['SNP.Name']==rmStartSNP].index[0]
				rmEndLoc = subDF[subDF['SNP.Name']==rmEndSNP].index[0]
				if workingLine>0 and rmMkrDF.loc[workingLine-1,'endLoc']>=rmStartLoc:
					rmMkrDF.loc[workingLine-1,'endLoc'] = rmEndLoc
					rmMkrDF.loc[workingLine-1,'endSNP'] = rmEndSNP
				else:
					rmMkrDF.loc[workingLine] = rmStartLoc,rmEndLoc,rmStartSNP,rmEndSNP
					workingLine +=1
			rmMkrsList = []
			for idx in rmMkrDF.index:
				rmMkrsList.extend(list(range(int(rmMkrDF.loc[idx,'startLoc']),int(rmMkrDF.loc[idx,'endLoc'])+1)))
			for k,g in groupby(enumerate(set(subDF.index)-set(rmMkrsList)), lambda x:x[0]-x[1]):
				includeMkrsList.extend(list(map(itemgetter(1),g)))
	length = len(includeMkrsList)
	if length<=int(maxNumsOfMkrs):
		selectGT = sampleDF.ix[includeMkrsList,]
	else:
		selectIndex = random.sample(includeMkrsList,maxNumsOfMkrs)
		selectGT = sampleDF.ix[selectIndex,]
	AA_baf = selectGT.ix[(selectGT['G1']=='A') & (selectGT['G2']=='A'),'BAF']
	BB_baf = selectGT.ix[(selectGT['G1']=='B') & (selectGT['G2']=='B'),'BAF']
	AB_baf = selectGT.ix[(selectGT['G1']=='A') & (selectGT['G2']=='B'),'BAF']
	
	delta_aa = np.percentile(AA_baf,float(perc))
	delta_bb = np.percentile(BB_baf,float(perc))
	sigma_ab = np.std(AB_baf)
	return([delta_aa,delta_bb,sigma_ab])

#######################################
## simulation
#######################################
def simData(seedNum,pediFilename,cleanMarkerFile,parentalHapDir,pennCNVDir,sampleList,dataSource,outputFilename,simLen,delta_aa,delta_bb,sigma_ab,pe1,simSDs,numOfSims):
	simDict = {}
	numpy.random.seed(seed=seedNum)
	## get input list information:
	## readin the lengthList file for the number of events and the length of each event.
	## list all files under the inputFilePath directory.
	allFiles = [f.strip() for f in open(sampleList).readlines()]
	## readin the parameters file:
	#delta_aa,delta_bb,sigma_ab,pe1 = map(float,modelParams.split(','))
	param = [float(delta_aa),float(delta_bb),float(sigma_ab)]
	pe1 = float(pe1)
	sd1,sd2 = map(float,simSDs.split('-'))
	###################################################################
	## get trioDict form myLib.py
	trioDict = getTrioDict(pediFilename, dataSource)
	markersDF = pd.read_table(cleanMarkerFile)
	## generate DataFrame
	## select one file at a time
	for iter in range(int(numOfSims)):
		while True:
			kidName = np.random.choice(list(allFiles))
			chrId = str(np.random.choice(range(1,23)))
			parentsNames_pair = trioDict[kidName][:2] # this generates a tuple
			pennCNVFilename = pennCNVDir + kidName + '.adjusted.rawcnv'
			rmMkrs = getParentalBadRegions(pennCNVFilename,parentsNames_pair)
			includeMkrsList,subDF = getNormalChrRegions(markersDF,chrId,rmMkrs)
			dupInfo = getNormalRegionOutput(includeMkrsList,subDF,int(simLen))
			if dupInfo:
				break
		#['rs2111400', 45522517, 'rs739012', 45850274]
		dadName,momName = parentsNames_pair
		startSNP = dupInfo[0]
		endSNP = dupInfo[2]
		eventLength = dupInfo[3]-dupInfo[1]+1
		#######################################################################
		## simulate BAF for offsprings
		dadFilename = parentalHapDir+'chr'+chrId+'/'+dadName+'_chr'+chrId+'.txt'
		momFilename = parentalHapDir+'chr'+chrId+'/'+momName+'_chr'+chrId+'.txt'
		dadPhased_df = pd.read_table(dadFilename,index_col=0,names=['SNP.Name','g1','g2'])
		momPhased_df = pd.read_table(momFilename,index_col=0,names=['SNP.Name','g1','g2'])
		
		f1 = dadPhased_df['g1']
		f2 = dadPhased_df['g2']
		m1 = momPhased_df['g1']
		m2 = momPhased_df['g2']
		hapDF = pd.concat([f1,f2,m1,m2],axis=1)
		hapDF.columns=['f1','f2','m1','m2']
		selected_hapDF = hapDF.loc[startSNP:endSNP,:]
		## check parental haplotypes homozygousity.
		f_homo_rate = sum(selected_hapDF['f1']==selected_hapDF['f2'])/float(simLen)
		m_homo_rate = sum(selected_hapDF['m1']==selected_hapDF['m2'])/float(simLen)

		fi = random.choice(['f1','f2']) # 'f1' or 'f2'
		mi = random.choice(['m1','m2']) # 'm1' or 'm2'
		# randomly pick one haplotype from each parent 
		h3 = random.choice(['f1','f2','m1','m2'])
		s = h3[0]
		if s=='f':
			if h3==fi:
				u = 0
			elif h3!=fi:
				u = 1
		elif s=='m':
			if h3==mi:
				u = 0
			elif h3!=mi:
				u = 1

		dupDF = selected_hapDF[[fi,mi,h3]]
		dupDF = dupDF.replace('-', np.nan)
		dupDF = dupDF.dropna()
		gtMap = {'AAA':'AAA','AAB':'AAB','ABA':'AAB','BAA':'AAB', 'ABB':'ABB','BAB':'ABB','BBA':'ABB','BBB':'BBB'}
		dupDF['gt'] = dupDF.apply(lambda x: gtMap[''.join(x)], axis=1)
		
		######################################################################
		## simulation:
		kidGT = simulate_gt(dupDF,sd1,sd2)
		#############################################################################################
		## calculate the pu1 and record if it is correct: estimate u --> p
		dat = pd.concat([Series(kidGT,index=dupDF.index,name='baf'), selected_hapDF],axis=1)
		p = dupCalc.calPu1_direct(dat, param, pe1)
		# f,m,h3,s,u,...
		outputList = [fi,mi,h3, s,str(u), str(p),kidName,dadName,momName, startSNP,endSNP, str(eventLength), simLen, str(f_homo_rate), str(m_homo_rate)]
		simDict[iter] = outputList
	print(len(simDict))
	simDF = pd.DataFrame.from_dict(simDict, orient='index')
	simDF.columns=['fi','mi','h3', 's','u', 'estimated_u','kidName','dadName','momName', 'startSNP','endSNP', 'eventLength', 'numOfMarkers', 'homoRate_father', 'homoRate_mother']
	simDF[['u','estimated_u','eventLength','numOfMarkers','homoRate_father','homoRate_mother']] = simDF[['u','estimated_u','eventLength','numOfMarkers','homoRate_father','homoRate_mother']].apply(pd.to_numeric)
	return(simDF)

#############################################################
###### calculate misMatch rate
###### translate from calMismatchRate.r
############################################################
def calMisCalRate(df,uEstThres,homoRateThres):
	df['diff'] = abs(df['u']-df['estimated_u'])
	df['homo_parent'] = np.where(df['s']=='f', df['homoRate_father'],df['homoRate_mother'])
	misMatchDF = df[(df['diff']>uEstThres) & (df['homo_parent']<homoRateThres)]
	r = (misMatchDF.shape[0])/(df.shape[0])
	return(r)


#######################################
## simulation, generate random seed as time
#######################################
def simData_timeRandom(pediFilename,cleanMarkerFile,parentalHapDir,pennCNVDir,sampleList,dataSource,outputFilename,simLen,delta_aa,delta_bb,sigma_ab,pe1,simSDs,numOfSims):
	simDict = {}
	random.seed(datetime.now())
	## get input list information:
	## readin the lengthList file for the number of events and the length of each event.
	## list all files under the inputFilePath directory.
	allFiles = [f.strip() for f in open(sampleList).readlines()]
	## readin the parameters file:
	#delta_aa,delta_bb,sigma_ab,pe1 = map(float,modelParams.split(','))
	param = [float(delta_aa),float(delta_bb),float(sigma_ab)]
	pe1 = float(pe1)
	sd1,sd2 = map(float,simSDs.split('-'))
	###################################################################
	## get trioDict form myLib.py
	trioDict = getTrioDict(pediFilename, dataSource)
	markersDF = pd.read_table(cleanMarkerFile)
	## generate DataFrame
	## select one file at a time
	for iter in range(int(numOfSims)):
		while True:
			kidName = np.random.choice(list(allFiles))
			chrId = str(np.random.choice(range(1,23)))
			parentsNames_pair = trioDict[kidName][:2] # this generates a tuple
			pennCNVFilename = pennCNVDir + kidName + '.adjusted.rawcnv'
			rmMkrs = getParentalBadRegions(pennCNVFilename,parentsNames_pair)
			includeMkrsList,subDF = getNormalChrRegions(markersDF,chrId,rmMkrs)
			dupInfo = getNormalRegionOutput(includeMkrsList,subDF,int(simLen))
			if dupInfo:
				break
		#['rs2111400', 45522517, 'rs739012', 45850274]
		dadName,momName = parentsNames_pair
		startSNP = dupInfo[0]
		endSNP = dupInfo[2]
		eventLength = dupInfo[3]-dupInfo[1]+1
		#######################################################################
		## simulate BAF for offsprings
		dadFilename = parentalHapDir+'chr'+chrId+'/'+dadName+'_chr'+chrId+'.txt'
		momFilename = parentalHapDir+'chr'+chrId+'/'+momName+'_chr'+chrId+'.txt'
		dadPhased_df = pd.read_table(dadFilename,index_col=0,names=['SNP.Name','g1','g2'])
		momPhased_df = pd.read_table(momFilename,index_col=0,names=['SNP.Name','g1','g2'])
		
		f1 = dadPhased_df['g1']
		f2 = dadPhased_df['g2']
		m1 = momPhased_df['g1']
		m2 = momPhased_df['g2']
		hapDF = pd.concat([f1,f2,m1,m2],axis=1)
		hapDF.columns=['f1','f2','m1','m2']
		selected_hapDF = hapDF.loc[startSNP:endSNP,:]
		## check parental haplotypes homozygousity.
		f_homo_rate = sum(selected_hapDF['f1']==selected_hapDF['f2'])/float(simLen)
		m_homo_rate = sum(selected_hapDF['m1']==selected_hapDF['m2'])/float(simLen)
			
		fi = random.choice(['f1','f2']) # 'f1' or 'f2'
		mi = random.choice(['m1','m2']) # 'm1' or 'm2'
		# randomly pick one haplotype from each parent 
		h3 = random.choice(['f1','f2','m1','m2'])
		s = h3[0]
		if s=='f':
			if h3==fi:
				u = 0
			elif h3!=fi:
				u = 1
		elif s=='m':
			if h3==mi:
				u = 0
			elif h3!=mi:
				u = 1
		
		dupDF = selected_hapDF[[fi,mi,h3]]
		dupDF = dupDF.replace('-', np.nan)
		dupDF = dupDF.dropna()
		gtMap = {'AAA':'AAA','AAB':'AAB','ABA':'AAB','BAA':'AAB', 'ABB':'ABB','BAB':'ABB','BBA':'ABB','BBB':'BBB'}
		dupDF['gt'] = dupDF.apply(lambda x: gtMap[''.join(x)], axis=1)

		######################################################################
		## simulation:
		kidGT = simulate_gt(dupDF,sd1,sd2)
		#############################################################################################
		## calculate the pu1 and record if it is correct: estimate u --> p
		dat = pd.concat([Series(kidGT,index=dupDF.index,name='baf'), selected_hapDF],axis=1)
		p = dupCalc.calPu1_direct(dat, param, pe1)
		# f,m,h3,s,u,...
		outputList = [fi,mi,h3, s,str(u), str(p),kidName,dadName,momName, startSNP,endSNP, str(eventLength), simLen, str(f_homo_rate), str(m_homo_rate)]
		simDict[iter] = outputList
	print(len(simDict))
	simDF = pd.DataFrame.from_dict(simDict, orient='index')
	simDF.columns=['fi','mi','h3', 's','u', 'estimated_u','kidName','dadName','momName', 'startSNP','endSNP', 'eventLength', 'numOfMarkers', 'homoRate_father', 'homoRate_mother']
	simDF[['u','estimated_u','eventLength','numOfMarkers','homoRate_father','homoRate_mother']] = simDF[['u','estimated_u','eventLength','numOfMarkers','homoRate_father','homoRate_mother']].apply(pd.to_numeric)
	return(simDF)

###############################################
## using unphased parental data
###############################################
#1. extract duplication region genotypes (offspring,father,mother)
# define offspring genotype (in the duplication region)
def getGenotypes(dupData_kid,dupData_dad,dupData_mom):
	kid_tmp = {}
	dad_tmp = {}
	mom_tmp = {}
	for i in dupData_kid.index:
		g1 = dupData_kid.loc[i,'G1']
		g2 = dupData_kid.loc[i,'G2']
		b = dupData_kid.loc[i,'BAF']
		if g1=="-" and g2=="-":
			dupG = "--"
		elif g1==g2:
			if b<=0.5:
				dupG = "AAA"
			elif b>0.5:
				dupG = "BBB"
		else:
			if b<=0.5:
				dupG = "AAB"
			elif b>0.5:
				dupG = "ABB"
		kid_tmp[i] = dupG
		pf_G = dupData_dad.loc[i,'G1']+dupData_dad.loc[i,'G2']
		pm_G = dupData_mom.loc[i,'G1']+dupData_mom.loc[i,'G2']
		dad_tmp[i] = pf_G
		mom_tmp[i] = pm_G

	kid_Gdata = Series(kid_tmp)
	dad_Gdata = Series(dad_tmp)
	mom_Gdata = Series(mom_tmp)
	
	dat = pd.concat([kid_Gdata,dad_Gdata,mom_Gdata], axis=1)
	dat.columns = ['kid','dad','mom']
	dat = dat.replace(to_replace="--",value=np.nan)
	# delete any SNP if there is a missing value:
	dat = dat.dropna(axis=0,how='any')
	L = len(dat)
	return(dat,L)

#2. the voting algorithm:
# create voting counters
def votingCount(dat):
	u0sF = 0; u1sF = 0;
	u0sM = 0; u1sM = 0;
	for snp in dat.index:
		kidG = dat.loc[snp,'kid']
		dadG = dat.loc[snp,'dad']
		momG = dat.loc[snp,'mom']
		if dadG=='AA' and momG=='AA':
			if kidG=='AAA':
				u0sF+=1; u1sF+=1; u0sM+=1; u1sM+=1
		elif dadG=='AA' and momG=='AB':
			if kidG=='AAA':
				u0sF+=1; u1sF+=1; u0sM+=1
			elif kidG=='AAB':
				u0sF+=1; u1sF+=1; u1sM+=1
			elif kidG=='ABB':
				u0sM+=1
		elif dadG=='AA' and momG=='BB':
			if kidG=='AAB':
				u0sF+=1; u1sF+=1
			elif kidG=='ABB':
				u0sM+=1; u1sM+=1
		elif dadG=='AB' and momG=='AA':
			if kidG=='AAA':
				u0sF+=1; u0sM+=1; u1sM+=1
			elif kidG=='AAB':
				u1sF+=1; u0sM+=1; u1sM+=1
			elif kidG=='ABB':
				u0sF+=1
		elif dadG=='AB' and momG=='AB':
			if kidG=='AAA':
				u0sF+=1; u0sM+=1
			elif kidG=='AAB':
				u0sF+=1; u1sF+=1; u0sM+=1; u1sM+=1
			elif kidG=='ABB':
				u0sF+=1; u1sF+=1; u0sM+=1; u1sM+=1
			elif kidG=='BBB':
				u0sF+=1; u0sM+=1
		elif dadG=='AB' and momG=='BB':
			if kidG=='AAB':
				u0sF+=1
			elif kidG=='ABB':
				u1sF+=1; u0sM+=1; u1sM+=1
			elif kidG=='BBB':
				u0sF+=1; u0sM+=1; u1sM+=1
		elif dadG=='BB'and momG=='AA':
			if kidG=='AAB':
				u0sM+=1; u1sM+=1
			elif kidG=='ABB':
				u0sF+=1; u1sF+=1
		elif dadG=='BB' and momG=='AB':
			if kidG=='AAB':
				u0sM+=1
			elif kidG=='ABB':
				u0sF+=1; u1sF+=1; u1sM+=1
			elif kidG=='BBB':
				u0sF+=1; u1sF+=1; u0sM+=1
		elif dadG=='BB' and momG=='BB':
			if kidG=='BBB':
				u0sF+=1; u1sF+=1; u0sM+=1; u1sM+=1
	return(u0sF,u1sF,u0sM,u1sM)

############# calculating PHASE CONCORDANCE for a region with a sample.
def getPhaseConcordanceScore(hapLOH_Dir,sampleID,markerFilename,startSNP,endSNP):
        count = 1
        markerDict = {}
        for line in open(markerFilename).readlines():
                SNPname = line.strip()
                markerDict[SNPname] = count
                count+=1
        startID = markerDict[startSNP]
        endID = markerDict[endSNP]

        informativeMarkerFilename = hapLOH_Dir + sampleID + '.informative'
        informativeMarkerTmp = open(informativeMarkerFilename).readline().strip().rsplit()
        informativeMarkerIdx = [int(x) for x in informativeMarkerTmp]
        print(len(informativeMarkerIdx))

        locatorStart = bisect.bisect(informativeMarkerIdx,startID)
        locatorEnd = bisect.bisect(informativeMarkerIdx,endID)
        print(locatorStart,locatorEnd)

        concordanceFilename = hapLOH_Dir + sampleID + '.switch_enumeration'
        concordanceTmp = open(concordanceFilename).readline().strip().rsplit()
        concordanceIdx = [int(x) for x in concordanceTmp]
        print(len(concordanceIdx))

        focusRegion = concordanceIdx[locatorStart:(locatorEnd+1)]
        print(len(focusRegion))
        try:
                phaseConcordance = float(sum(focusRegion))/len(focusRegion)
        except ZeroDivisionError:
                phaseConcordance = 'NaN'
                print('check %s!' % sampleID)
        print(phaseConcordance)
        return phaseConcordance

####### calculating BAF (heterozygous markers only) standard deviation etc.
def calBAF_stat(filename,missingLabel):
        ## This is for checking the data quality:
        # calculating the mean, median, and variation of BAF values of heterozygous markers in a region

        # formatted samples from step1
        # different directory for blood samples and nonBlood samples
        #SNP.Name       Chr     Position        BAF     LRR     G1      G2
        #rs12255619     10      88481   0.0015  -0.0789 A       A
        #rs11591988     10      116070  0.9999  0.2305  G       G
        #rs4508132      10      121636  0.4545  -0.0168 -       -
        df = pd.read_table(filename,na_values=missingLabel)
        df1 = df.dropna()
        het_index = df1[ df1['G1']!=df1['G2'] ].index
        df2 = df1.loc[het_index]
        x = np.array(df2['BAF'])
        print('BAF_mean (het markers) = %f' % numpy.mean(x))
        print('BAF_median (het markers) = %f' % numpy.median(x))
        print('BAF_std (het markers) = %f' % numpy.std(x))
        return([np.mean(x),np.median(x),np.std(x)])

####### calculating LRR standard deviation etc.
def calLRR_stat(filename,missingLabel):
        df = pd.read_table(filename,na_values=missingLabel)
        df1 = df.dropna()
        y = np.array(df1['LRR'])
        print('LRR_mean = %f' % numpy.mean(y))
        print('LRR_median = %f' % numpy.median(y))
        print('LRR_std = %f' % numpy.std(y))
        return([np.mean(y),np.median(y),np.std(y)])

####### calculating missing genotype call rate:
## This is for checking the data quality:
# calculating the missing genotype rate for a single sample.
## The input data format requirement:
#SNP.Name       Chr     Position        BAF     LRR     G1      G2
#rs12255619     10      88481   0.0015  -0.0789 A       A
#rs11591988     10      116070  0.9999  0.2305  G       G
#rs4508132      10      121636  0.4545  -0.0168 -       -
def calMissingCallRate(filename, missingLabel):
        '''
        filename contains the sample data information following the format requirement.
        missingLabel is the one character need to be considered as missing in the two genotype columns G1 & G2.
        '''
        df = pd.read_table(filename, index_col=0,usecols=[0,5,6],na_values=missingLabel)
        totalCount = len(df.index)
        missingCount = max(df.isnull().sum())
        missingCallRate = missingCount/totalCount
        print('Sample %s has missing genotype call rate at %f' % (filename, missingCallRate))
        return(missingCallRate)
