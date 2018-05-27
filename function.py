#!/usr/bin/python3
# -*",-coding:Latin-1 -utf8",
import os
import re
import sys
from decimal import *
import FastaSeq
import random
import math
import collections 
import inFile_related_function
import blast_related


"""
	This module retrieves all general functions use in main program
"""
def RevCom(sequence):
	revSeq = []
	compRevSeq = []
	seq = ""
	nucleotides = list(sequence)
	for curent_nucleotide in nucleotides:
		if curent_nucleotide =="T":
			revSeq.append("A")
		if curent_nucleotide =="A":
			revSeq.append("T")
		if curent_nucleotide =="G":
			revSeq.append("C")
		if curent_nucleotide =="C":
			revSeq.append("G")
	revSeq.reverse()
	for finalRevSeq in revSeq:
		seq+= finalRevSeq
	return seq

 
def setExonList (exonFind,filename="none"):
	fileOUt= ".\\result\\"+filename+"_"+"exon_list.txt"	
	result = open(fileOUt,"w")
	result.write("**** Isoform:")
	filename = filename.replace("fa","")
	result.write(filename)
	result.write(" ****")
	result.write("\n")
	for index,seq in exonFind.items():
		result.write(">exon ")
		result.write(str(index))
		result.write("\n")
		result.write(str(seq))
		result.write("\n")
	result.close()
	return fileOUt

def CountNucleotide(sequence):
	#intialization of counter
	countA,countT,countC,countG = 0,0,0,0
	nucleotides = list(sequence)
	for curent_nucleotide in nucleotides:
		if curent_nucleotide =="A":
			countA = countA+1
		if curent_nucleotide =="T":
			countT = countT+1
		if curent_nucleotide =="C":
			countC = countC+1
		if curent_nucleotide =="G":
			countG = countG+1
	return countA,countT,countC,countG


def setCircosKaryo(blastRes,filename="none"):
	foldername = inFile_related_function.createResultFolder("karyotype")
	fileOUt= foldername+"\\"+filename+"_"+"karyotype.exon.txt"
	result = open(fileOUt,"w")
	result.write("#chr - ID LABEL START END COLOR")
	result.write("\n")
	for name,arrayOfRes in blastRes.items():
		for index in range (len(arrayOfRes)):
			tempArray = arrayOfRes[index].split('\t')
			qStart = int(tempArray[6])
			qStop =	int(tempArray[7])		
			result.write("chr - axis")
			result.write(str(index))
			result.write(" exon")
			result.write(str(index))
			result.write(" ")
			result.write(str(qStart))
			result.write(" ")
			result.write(str(qStop))	
			if index %2 == 0:
				result.write(" blue")
			else :
				result.write(" green")
			result.write("\n")
	result.close()
	return fileOUt

"""
	gc content 
"""
def GCcontent(sequence):
	totalNuc = 0
	GCcount = 0
	nucleotides = list(sequence)
	for curent_nucleotide in nucleotides:
		if curent_nucleotide =="C":
			GCcount = GCcount+1
		if curent_nucleotide =="G":
			GCcount = GCcount+1
	totalNuc = len(sequence)
	GCcount= Decimal(GCcount)
	totalNuc= Decimal(totalNuc)
	GCrate = 100*GCcount/totalNuc
	GCrate.quantize(Decimal('.00001'), rounding=ROUND_HALF_UP)
	GCrate = round(GCrate,6)
	return GCrate

"""
	evaluation of tmRate 

"""
def TmCheck(sequence):
	countA,countT,countC,countG = 0,0,0,0
	countA,countT,countC,countG = CountNucleotide(sequence)
	TmRate = ((4*(countG+countC))+ ((2*(countA+countT))-4))
	return TmRate




"""
	TODO check the 3' complementary  

"""


def primerCdna(exonFind):
	firstSeq = None
	forward = ""
	reverse = ""
	for index,seq in exonFind.items():		
			forward,reverse = FastaSeq.getPrimers(seq,cuttOff=20)
			print (seq,"\n",forward,reverse)
			# firstSeq = seq



def setPrimerOutput(finalArray,filename="none"):
	fileOUt = ".\\result\\"+filename+"_"+"primer_list.txt"	
	result = open(fileOUt,"w")	
	result.write("**** Isoform ")
	filename = filename.replace("fa","")
	result.write(filename)
	result.write(" ****")
	result.write("\n")
	result.write("#exon number, primerFW, primerRV")
	result.write("\n")
	for name,seq2input in finalArray.items():
		amorceFW,amorceRV = FastaSeq.getPrimers(seq2input)		
		result.write(str(name))	
		# result.write(" ")	
		# result.write(seq2input)
		result.write(",")		
		if str(amorceFW).startswith("(\"not available ") or str(amorceFW).startswith("(\'not available ") :
			amorceFW = "not available"		
		if str(amorceRV).startswith("(\"not available ") or str(amorceRV).startswith("(\'not available ")  :
			amorceRV = "not available"
		result.write(str(amorceFW))			
		result.write(",")	
		result.write(str(amorceRV))
		result.write("\n")
	result.close()
	return fileOUt




def findExonSetPrimer(genomicFile,cdnaFile,blastRes):
	genomicDict =  inFile_related_function.ExtractFasta(genomicFile)
	cdnaDict =  inFile_related_function.ExtractFasta(cdnaFile)
	finalArray = {}
	outExon = ""
	for genomicName,genomicSeq in genomicDict.items():
		print ("In progress ...\n")
		for cdnaName,cdnaSeq in cdnaDict.items():
			if  len(cdnaSeq) > len(genomicSeq):
				print("You have input the wrong genomic file, cdna sequence is longer than genomic sequence.\n Please check your input files and try again")
				exit ()
			finalArray = setExonDic(blastRes,cdnaSeq)
			 
	outfile = os.path.basename(cdnaFile) # To keep the name of the file
	outExon = setExonList (exonFind,outfile)
	outFile = setPrimerOutput(finalArray,outfile)
	print ("These files\n","*",outExon,"\n*",outFile,"\nhave been created.\nThey contain the exon list and the list of finding amorces.")






def setExonDic(blastRes,cdnaSeq):
	exonFind = {}
	startPos = []
	temp = []
	startStop = {}
	for name,arrayOfRes in blastRes.items():
		for index in range (len(arrayOfRes)):
			tempArray = arrayOfRes[index].split('\t')
			qStart = int(tempArray[6])
			qStop =	int(tempArray[7])
			temp.append(qStart) # I select all start position						
			startPos = sorted(temp) 
			startStop[qStart]= qStop 
			i=1
		for start in startPos:
			if start == 1 :
				sequence = cdnaSeq[0:startStop.get(start)+10]
				exonFind[i] = sequence
			elif start == startPos[-1]:
				sequence = cdnaSeq[start:startStop.get(start)]
				exonFind[i] = sequence
			else :
				sequence = cdnaSeq[start:startStop.get(start)+10]
				exonFind[i] = sequence
			i+=1		
	return exonFind


def runMain(genomicFile,cdnaFile):
	cdnaDict =  inFile_related_function.ExtractFasta(cdnaFile)
	genomicDict =  inFile_related_function.ExtractFasta(genomicFile)

# -- Creation of output files to fit in circos program
#   Thanks to blast alignment program a list of exon is extracted from genomic and cdna files 

		# -- Run blast 

	blastResultfile = blast_related.runBlastN(cdnaFile,genomicFile)
	blastRes = blast_related.parseBlast(blastResultfile)



	
	for genomicName,genomicSeq in genomicDict.items():
		print ("In progress ...\n")
		for cdnaName,cdnaSeq in cdnaDict.items():
			if  len(cdnaSeq) > len(genomicSeq): # check of input files
				print("You have input the wrong genomic file, cdna sequence is longer than genomic sequence.\n Please check your input files and try again")
				exit ()
		# -- Extract exon from blast result	

			finalArray = setExonDic(blastRes,cdnaSeq) 
		# -- creation of amorce list extract from exon values 

	
# -- creation of outfile 

	outfile = os.path.basename(cdnaFile) # To keep the name of the file
	outExon = setExonList (finalArray,outfile)
	setPrimerOutput(finalArray,outfile)
	outCircos = setCircosKaryo(blastRes,outfile)

	print ("This file\n","*",outExon,"\nhas been created.\nIt contains the exon list.")
	return outExon

def endOfprogam(quit=None):
	while quit == None or quit == "R" or quit == "rerun" or quit == "r" or quit == "RERUN":
		quit = input("Press Q to exit this program or R or rerun this program again:\n")
		if quit == "R" or quit == "rerun" or quit == "r" or quit == "RERUN":
			genomicFile = input("Enter the genomic file name as it is written in the current folder.\n")
			cdnaFile = input("Enter the cdna file name as it is written in the current folder.\n")
			runMain(genomicFile,cdnaFile)
		elif quit == "Q" or quit == "q": 
			print("End of the program, have a nice day :D.")
			sys.exit()
		else:
			print(quit,"is not a correct option! End of the program have a nice day!\n")
			sys.exit()
