#!/usr/bin/python3.4
# -*",-coding:Latin-1 -utf8",
import os
import re
import sys
import FastaSeq
import shutil



# -- set up final exon list 
# It pastes together all exon list files
def finalExonList(filename="Test",*exonListFile):
	fileOUt = ".\\result\\"+filename+"_finalExonList.txt"
	exonReport = open(fileOUt,"w")
	for exonFile in exonListFile:
		exonFile = ".\\result\\"+exonFile
		shutil.copyfileobj(open(exonFile, 'r'), exonReport)
	exonReport.close()
	return fileOUt

def primerInCommon(exonReport):
	Fwreport = {} # dictionnary with isoform name as key and Forward primer list as values
	Rvreport = {} # dictionnary with isoform as key and reverse primer list as values
	FwCount = {} # dictionnary with isoform as key and the number of this isoform in the input file
	RvCount = {} # dictionnary with isoform as key and the number of this isoform in the input file
	primerFw = []
	primerRv = []
	exonReport = exonReport
	name = ""
	try :
		fileIn = open(exonReport,"r")
	except FileNotFoundError:
		print("Your file",exonReport," does not exist !")
		exit()
	for exon in fileIn:
		exon =exon.strip()	
		# print(exon)
		if not exon.startswith('>exon') and not exon.startswith('***'):
			amorceFW,amorceRV = FastaSeq.getPrimers(exon)
			if str(amorceFW).startswith("(\"not available") or str(amorceFW).startswith("(\'not available"):
				amorceFW = 'not available'
			if str(amorceRV).startswith("(\"not available") or str(amorceRV).startswith("(\'not available"):
				amorceRV = 'not available'
			primerFw.append(amorceFW)
			primerRv.append(amorceRV)
		if exon.startswith('***'):
			exon = exon.replace("**** Isoform:","")
			exon = exon.replace(" ****","")
			name = exon
			Fwreport[name] = primerFw
			Rvreport[name] = primerRv
	fileIn.close()
	# -- count the number of primers and put it in a dictionnary 
	FwCount = countPrimers(primerFw) 
	RvCount = countPrimers(primerRv)
	return Fwreport, Rvreport, FwCount, RvCount
	
def countPrimers(arrayOfPrimer):
    countPrimer = {}.fromkeys(set(arrayOfPrimer),0)
    for primer in arrayOfPrimer:
        countPrimer[primer] += 1
    return countPrimer



file1 = "TraesCS2A01G514800.1.fa_exon_list.txt"
file2 = "TraesCS2A01G514800.2.fa_exon_list.txt"
file3 = "TraesCS2B01G543100.1.fa_exon_list.txt"
file4 = "TraesCS2D01G516300.1.fa_exon_list.txt"
file5 = "TraesCS2D01G516300.2.fa_exon_list.txt"


exonReport = finalExonList("EME1A_AT2G21800.2",file1,file2,file3,file4,file5)
primerInCommon(exonReport)


def setCsvReport(reportFw,countFw,reportRv,countRv,filename="Test"):
	fileOUt = ".\\result\\"+filename+"_finalReport.txt"
	PrimerRate,primer,numOfIsoform = 0,0,0
	numOfIsoform = len(reportFw.keys())
	csvReport = open(fileOUt,"w")
	csvReport.write("###forward , reverse, in common, not in common(isoform)")
	csvReport.write("\n")
	for name,arrayOfprimer in reportFw.items():
		for seq in arrayOfprimer:
			if not seq.startswith("not available"):
				PrimerRate = countFw.get(seq)	 			
				if PrimerRate == numOfIsoform :
					csvReport.write(seq)
					csvReport.write(",N\\A,yes,all")
					csvReport.write("\n")
					# print(seq,",N\\A,yes,all")
				else:
					csvReport.write(seq)
					csvReport.write(",N\\A,no,")
					csvReport.write(name)
					csvReport.write("\n")
					# print(seq,",N\\A,no,",name)
	for name,arrayOfprimer in reportRv.items():
		for seq in arrayOfprimer:
			if not seq.startswith("not available"):
				PrimerRate = countRv.get(seq)	 			
				if PrimerRate == numOfIsoform :
					csvReport.write("N\\A,")
					csvReport.write(seq)
					csvReport.write(",yes,all")
					csvReport.write("\n")
					# print("N\\A,",seq,",yes,all")
				else:
					csvReport.write("N\\A,")
					csvReport.write(seq)
					csvReport.write(",no,")
					csvReport.write(name)
					csvReport.write("\n")

	csvReport.close()
	return fileOUt
		

Fwreport, Rvreport, FwCount, RvCount = primerInCommon(exonReport)

setCsvReport(Fwreport,FwCount,Rvreport,RvCount,"EME1A_AT2G21800.2")