#!/usr/bin/python3
# -*",-coding:Latin-1 -utf8",
import os
import re
import sys
from decimal import *
import FastaSeq
import function
import inFile_related_function
import output_interface_related_function


"""
This module is the main module

"""



#-- Initialization of variables 

finalArray = {}
outExon = ""

# -- Creation of result files folder 

inFile_related_function.createResultFolder()

# -- output_interface_related_function.firstPage ("Hello!\nPlease enter the genomic file name as it is written.\n Please pay attention that the input file and this program must be in the same folder.\n")# This to output the message on text box

# -- Extract of genomic and cdna sequences from input and creation of dictionnary with sequences as value and name of sequences as keys

genomicFile = input("Hello!\nPlease enter the genomic file name as it is written.\n Please pay attention that the input file and this program must be in the same folder.\n") 

genomicDict =  inFile_related_function.ExtractFasta(genomicFile)

cdnaFile = input("Please enter the cdna file name as it is written.\n Pay attention that the input file and this program must be in the same folder.\n")

cdnaDict =  inFile_related_function.ExtractFasta(cdnaFile)
"""
Creation of output to fit in circos program
a list of exon from input genomic and cdna 
creation of amorce list to fit to serial cloner program 

"""
for genomicName,genomicSeq in genomicDict.items():
	print ("In progress ...\n")
	for cdnaName,cdnaSeq in cdnaDict.items():
		if  len(cdnaSeq) > len(genomicSeq):
			print("You have input the wrong genomic file, cdna sequence is longer than genomic sequence.\n Please check your input files and try again")
			exit ()
		exonFind = function.getExon(genomicSeq,cdnaSeq)
		finalArray = function.setSequence(exonFind)

"""
set up of output files 

"""
outfile = os.path.basename(cdnaFile) # To keep the name of the file

outExon = function.setExonList (exonFind,outfile)
outAmorce = function.setAmorceOutput(finalArray,outfile)
outCircos = function.setCircosKaryo(exonFind,outfile)
print ("These files\n","*",outExon,"\n*",outAmorce,"\nhave been created.\nThey contain the exon list and the list of finding amorces.")

""" 
creation of new amorces list from output files and setup the end of the program

"""

amoreFromExonFile = input("Would you like to design, your amorce from output exon list (Yes/No)?\n")
if amoreFromExonFile == "Y" or amoreFromExonFile == "Yes" or amoreFromExonFile == "y" or amoreFromExonFile == "yes":
	exonSelected = ""
	numOfExon = input("Write the number of your selected exon:\n")
	exonSelected = "exon "+numOfExon
	function.amorceFromInputExon(exonSelected,outfile,cuttOff = 20)
else:
	function.endOfprogam()
