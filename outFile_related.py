#!/usr/bin/python3.4
# -*",-coding:Latin-1 -utf8",
import os
import re
import sys
import FastaSeq
import shutil
import glob



# -- group and set up a final amorce list from all input primer list 
#TODO add comment on each functions
def mergePrimers(folderName,primers,special=None):
	if special == None :
		fileOUt = ".\\"+folderName+"\\"+folderName+"_allPrimers.txt"
	elif special:
		fileOUt = ".\\"+folderName+"\\"+folderName+"_"+special+"Primers.txt"
	primerReport = open(fileOUt,"w")
	for primerList in primers:
		shutil.copyfileobj(open(primerList, 'r'), primerReport)
	primerReport.close()
	return fileOUt


def finalExonList(folderName,exonListFile):
	fileOUt = ".\\"+folderName+"\\"+folderName+"_allExons.txt"
	exonReport = open(fileOUt,"w")
	for exonFile in exonListFile:
		shutil.copyfileobj(open(exonFile, 'r'), exonReport)
	exonReport.close()
	return fileOUt

def primerInCommon(primerReport):
	Fwreport = {} # dictionnary with isoform name as key and Forward primer list as values
	exonFw = {} # 
	exonRv = {} # 
	Rvreport = {} # dictionnary with isoform as key and reverse primer list as values
	FwCount = {} # dictionnary with isoform as key and the number of this isoform in the input file
	RvCount = {} # dictionnary with isoform as key and the number of this isoform in the input file
	arrayFw = []
	arrayRv = []
	# exonReport = exonReport
	name = ""
	try :
		fileIn = open(primerReport,"r")
	except FileNotFoundError:
		print("Your file",primerReport," does not exist !")
		exit()
	for line in fileIn:
		line = line.strip()	
		
		if line.startswith('#') :
			pass
		elif line.startswith('**** '):
			line = line.replace("**** Isoform","")
			line = line.replace(" ****","")
			name =line	
			Fwreport[name] = arrayFw
			Rvreport[name] = arrayRv		
		else :
			temp = line.split(',')
			exon = temp[0]+name
			primerFw = temp[1]
			primerRv = temp[2]
			arrayFw.append(primerFw)
			arrayRv.append(primerRv)
			exonFw[primerFw] = exon
			exonRv[primerRv] = exon
	
	
	fileIn.close()
	# -- count the number of primers and put it in a dictionnary 
	FwCount = countPrimers(arrayFw) 
	RvCount = countPrimers(arrayRv)
	return Fwreport, Rvreport, FwCount, RvCount,exonFw,exonRv


def countPrimers(arrayOfPrimer):
    countPrimer = {}.fromkeys(set(arrayOfPrimer),0)
    for primer in arrayOfPrimer:
        countPrimer[primer] += 1
    return countPrimer


def setCsvReport(folderName,reportFw,countFw,reportRv,countRv,exonFw,exonRv,special=None):
	if special == None :
		fileOUt = ".\\"+folderName+"\\"+folderName+"_finalReport.txt"
	elif special:
		fileOUt = ".\\"+folderName+"\\"+folderName+"_"+special+"_finalReport.txt"
	
	PrimerRate,primer,numOfIsoform = 0,0,0
	numOfIsoform = len(reportFw.keys())
	csvReport = open(fileOUt,"w")
	csvReport.write("### exon isoform, forward, reverse, in common( yes or copy number)")
	csvReport.write("\n")
	for name,primerFw in reportFw.items():
		for seq in primerFw:			
			if not seq.startswith("not available"):
				exon = exonFw.get(seq)	
				csvReport.write(exon)
				csvReport.write(",")
				PrimerRate = countFw.get(seq)			 			
				if PrimerRate == numOfIsoform :
					csvReport.write(seq)
					csvReport.write(",N\\A,yes")
					csvReport.write("\n")
				else:
					csvReport.write(seq)
					csvReport.write(",N\\A,")
					csvReport.write(str(PrimerRate))
					csvReport.write("\n")
					# print(seq,",N\\A,no,",name)
	for name,primerRv in reportRv.items():
		for seq in primerRv:			
			if not seq.startswith("not available"):
				exon = exonRv.get(seq)	
				csvReport.write(exon)
				csvReport.write(",")
				PrimerRate = countRv.get(seq)	 			
				if PrimerRate == numOfIsoform :
					csvReport.write("N\\A,")
					csvReport.write(seq)
					csvReport.write(",yes")
					csvReport.write("\n")
				else:
					csvReport.write("N\\A,")
					csvReport.write(seq)
					csvReport.write(",")
					csvReport.write(str(PrimerRate))
					csvReport.write("\n")

	csvReport.close()
	return fileOUt

def generateReport(foldername) :
	file = glob.glob(".\\"+foldername+"\\*.fa_exon_list.txt")
	primers = glob.glob(".\\"+foldername+"\\*.fa_primer_list.txt")	
	genoPrimers = glob.glob(".\\"+foldername+"\\*.fa_genomic_primer_list.txt")
	exonReport = finalExonList(foldername,file)
	primerReport = mergePrimers(foldername,primers)
	genoPrimerReport = mergePrimers(foldername,genoPrimers,"genomic")
	Fwreport, Rvreport, FwCount, RvCount,exonFw,exonRv = primerInCommon(primerReport)
	genomicReport = setCsvReport(foldername,Fwreport,FwCount,Rvreport,RvCount,exonFw,exonRv)
	genoFwreport, genoRvreport, genoFwCount, genoRvCount,genoExonFw,genoExonRv = primerInCommon(genoPrimerReport)
	cdnaReport = setCsvReport(foldername,genoFwreport,genoFwCount,genoRvreport,genoRvCount,genoExonFw,genoExonRv,"genomic")
	return genomicReport,cdnaReport