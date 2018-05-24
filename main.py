#!/usr/bin/python3
# -*",-coding:Latin-1 -utf8",
import os
import re
import sys
from decimal import *
import FastaSeq
import function
import inFile_related_function
import blast_related


"""
This module is the main module

"""



#-- Initialization of variables, tuples and dictionnary

finalArray = {}
genomicDict = {}
cdnaDict = {}
genomicFile = ""
cdnaFile = ""
outExon = ""
outfile = ""
blastRes = ""
blastResultfile = ""

# -- Creation of result folder, it contains all output of this program 

inFile_related_function.createResultFolder()

# -- Extract of genomic and cdna sequences from input files. It creates a dictionnary with sequences as value and name of sequences as keys for genomic and cdna files. The name of genomic and cdna files are input from keyboard values 

genomicFile = input("Hello!\nPlease enter the genomic file name as it is written.\n Please pay attention that the input file and this program must be in the same folder.\n") 

genomicDict =  inFile_related_function.ExtractFasta(genomicFile)

cdnaFile = input("Please enter the cdna file name as it is written.\n Pay attention that the input file and this program must be in the same folder.\n")

cdnaDict =  inFile_related_function.ExtractFasta(cdnaFile)

# -- Creation of output files to fit in circos program
#   Thanks to blast alignment program a list of exon is extracted from genomic and cdna files 


		# -- Run blast 

blastResultfile = blast_related.runBlastN(cdnaFile,genomicFile)
blastRes = blast_related.parseBlast(blastResultfile)

# -- creation of amorce list extract from exon values 

for genomicName,genomicSeq in genomicDict.items():
	print ("In progress ...\n")
	for cdnaName,cdnaSeq in cdnaDict.items():
		if  len(cdnaSeq) > len(genomicSeq): # check of input files
			print("You have input the wrong genomic file, cdna sequence is longer than genomic sequence.\n Please check your input files and try again")
			exit ()
		# -- Extract exon from blast result	

		finalArray = function.setExonDic(blastRes,cdnaSeq) 
	
# -- creation of outfile 

outfile = os.path.basename(cdnaFile) # To keep the name of the file
outExon = function.setExonList (finalArray,outfile)
# outAmorce = function.setAmorceOutput(finalArray,outfile)
outCircos = function.setCircosKaryo(blastRes,outfile)

print ("This file\n","*",outExon,"\nhas been created.\nIt contains the exon list.")


# -- creation of new amorces list from output files 

amoreFromExonFile = input("Would you like to design, your amorce from output exon list (Yes/No)?\n")
if amoreFromExonFile == "Y" or amoreFromExonFile == "Yes" or amoreFromExonFile == "y" or amoreFromExonFile == "yes":
	exonSelected = ""
	numOfExon = input("Write the number of your selected exon:\n")
	exonSelected = "exon "+numOfExon
	function.amorceFromInputExon(exonSelected,outfile,cuttOff = 20)
else:
	function.endOfprogam()
