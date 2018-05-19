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
	print (exonFind)	
	fileOUt= ".\\result\\"+filename+"_"+"exon_list.txt"	
	result = open(fileOUt,"w")
	for index,seq in exonFind.items():
		name = index +1
		result.write("exon ")
		result.write(str(name))
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


def setCircosKaryo(exonFind,filename="none"):
	foldername = inFile_related_function.createResultFolder("karyotype")
	fileOUt= foldername+"\\"+filename+"_"+"karyotype.exon.txt"
	result = open(fileOUt,"w")
	result.write("#chr - ID LABEL START END COLOR")
	result.write("\n")
	for name,seq in exonFind.items():		
		result.write("chr - axis")
		result.write(str(name))
		result.write(" exon")
		result.write(str(name))
		result.write(" 0 ")
		size = 100*len (seq) # exon size is too short 
		result.write(str(size))
		
		if name % 2 == 0:
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
	# print ("le gc content est de",GCrate)
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
	check the 3' complementary  

"""

def setAmorceOutput(finalArray,filename="none"):
	fileOUt = ".\\result\\"+filename+"_"+"amorce_list.txt"
	result = open(fileOUt,"w")
	for name,seq2input in finalArray.items():
		amorceFW,amorceRV = FastaSeq.getPrimers(seq2input)
		result.write("Sequence ")
		result.write(str(name))
		result.write("\n")
		result.write(seq2input)
		result.write("\n")
		result.write("amorceFW ")
		result.write("\n")
		result.write(str(amorceFW))
		result.write("\n")
		result.write("amorceRV ")
		result.write("\n")
		result.write(str(amorceRV))
		result.write("\n")
	result.close()
	return fileOUt




def findExonSetAmorce (genomicFile,cdnaFile,blastRes):
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
	outFile = setAmorceOutput(finalArray,outfile)
	print ("These files\n","*",outExon,"\n*",outFile,"\nhave been created.\nThey contain the exon list and the list of finding amorces.")





def amorceFromInputExon(exonSelected,filename,cuttOff=20):
	exon2Add = ""
	sequence = ""
	add = None
	name = ""
	if inFile_related_function.isExist(exonSelected,filename):
		exon2Add = inFile_related_function.getExon(exonSelected,filename)
		sequence+= exon2Add
		name+= str(exonSelected)+"_"
	else:
		print("This exon does not exist\nEnd of the program")
		exit()
	while len(sequence) < 70 :
		numOfExon = input("the length of your selected exon is too short. Add another one.\n")
		exonSelected = "exon "+numOfExon
		if inFile_related_function.isExist(exonSelected,filename):
			exon2Add = inFile_related_function.getExon(exonSelected,filename)
			sequence+= exon2Add			
			name+= str(exonSelected)+"_"	
			# print(sequence)
		elif not inFile_related_function.isExist(exonSelected,filename):
			print("This exon does not exit\nEnd of the program")			
			exit()	
	if len(sequence) > 70:	
		add = input("Would you like to add an other exon? Yes/No.\n")	
		if add  == "Yes":
			numOfExon = input("Write the number of the exon that you want to add.\n")
			exonSelected = "exon "+numOfExon
			if inFile_related_function.isExist(exonSelected,filename):
				exon2Add = inFile_related_function.getExon(exonSelected,filename)
				sequence+= exon2Add			
				name+= str(exonSelected)+"_"
		elif add  == "No":
			fileOUt = ".\\result\\"+filename+"_"+"exonSelected_amorce.txt"
			result = open(fileOUt,"w")
			amorceFW,amorceRV = FastaSeq.getPrimers(sequence,cuttOff = 20)
			GCrate = GCcontent(amorceFW)
			print ("1:",GCrate)	
			GCrate = GCcontent(amorceRV)		
			print ("2:",GCrate)			
			result.write("Sequence ")
			result.write(name)
			result.write("\n")
			result.write(sequence)
			result.write("\n")
			result.write("amorceFw ")
			result.write("\n")
			result.write(amorceFW)
			result.write("\n")
			result.write("amorceRv")
			result.write("\n")
			result.write(amorceRV)
			result.write("\n")
			result.close()
			print("result file is created here:",fileOUt,"\n")
		else:
			Print("I don't understand your choice.")
			add  = input("Would you like to add an other exon? Yes/No.\n")






def endOfprogam(quit=None):
	while quit == None or quit == "R" or quit == "rerun" or quit == "r" or quit == "RERUN":
		quit = input("Press Q to exit this program or R or rerun this program again:\n")
		if quit == "R" or quit == "rerun" or quit == "r" or quit == "RERUN":
			genomicFile = input("Enter the genomic file name as it is written in the current folder.\n")
			cdnaFile = input("Enter the cdna file name as it is written in the current folder.\n")
			findExonSetAmorce (genomicFile,cdnaFile)
		elif quit == "Q" or quit == "q": 
			print("End of the program, have a nice day :D.")
			sys.exit()
		else:
			print(quit,"is not a correct option! End of the program have a nice day!\n")
			sys.exit()


def setExonDic(blastRes,cdnaSeq):
	exonFind = {}
	for name,arrayOfRes in blastRes.items():
		for index in range (len(arrayOfRes)):
			tempArray = arrayOfRes[index].split('\t')
			qStart = int(tempArray[6])
			qStop =	int(tempArray[7])
			sequence = cdnaSeq[qStart:qStop]			
			exonFind [index]= sequence 
	return exonFind
