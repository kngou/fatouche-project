#!/usr/bin/python3
# -*",-coding:Latin-1 -utf8",
import os
import re
import sys
from decimal import *
import FastaSeq
import function
import readCommandLine
import inFile_related_function
import blast_related
import glob
import getopt
import outFile_related
# from optparse import OptionParser


"""
This module is the main module

"""
# -- read input command line

folderName,master,genomicList,cdnaList,verbose = readCommandLine.commandLine()
# print(folderName,master,genomicFile,cdnaFile,verbose )


# sys.exit()
# -- Creation of result folder by default the folder name is "result" but it can be changed thanks to the command line
#    it will contain all output files of this program but circos related files 

inFile_related_function.createResultFolder(folderName)

"""
# -- Extraction of genomic and cdna sequences from input files. 
		It creates a dictionnary with sequences as value and name of sequences as keys for all input files.
			two options are available : 
				* master option : thanks to -m or -- master option in command line : 
					it searches all files with .fasta extension and put it on a list.
					Then it run the program on it. It assumes that the corresponding cdna file has the same name but with the .fa extension.
				*genomicName and cdnaName retrieve thanks to -gF -cF or cdnaFile and genomicFile

"""
if (master and  cdnaList) or (master and genomicList):
	print ("You can not choose these options together, please press -h to know the correct option.")
	sys.exit()
elif master:
	file = glob.glob('*.fasta')
	NumOfFile = len(file)
	allFileIn = 0
	for file in file :
		genomicFile = file
		cdnaFile = file.replace("fasta","fa")	
		function.runMain(genomicFile,cdnaFile,folderName) # output of exonlist from genomic and cdna files.
		allFileIn+=1

	if NumOfFile == allFileIn:
		genomicReport,cdnaReport = outFile_related.generateReport(folderName) # outputs final report with amorce list
elif genomicList and cdnaList:	
	allFileIn = 0
	genomicFile = str(genomicList).split(",")
	cdnaFile = str(cdnaList).split(",")
	genoSize = len(genomicFile)
	cdnaSize = len(cdnaFile)	
	if genoSize != cdnaSize:
		print("something went wrong with input files!")
		sys.exit()
	elif genoSize == cdnaSize:
		for index in range(genoSize):
			genomic = genomicFile[index].replace("[\'","")
			genomic= genomic.replace('\']','')
			cdna = cdnaFile[index].replace('[\'','')
			cdna = cdna.replace('\']','')
			# print(cdna,genomic,folderName) # output of exonlist from genomic and cdna files.
			function.runMain(genomic,cdna,folderName) # output of exonlist from genomic and cdna files.
		
		outfiles = glob.glob(".\\"+folderName+"\\*.fa_exon_list.txt")
		if len (outfiles) == genoSize:
			genomicReport,cdnaReport = outFile_related.generateReport(folderName) # outputs final report with amorce list
		elif genoSize == 1:
			genomicReport,cdnaReport = outFile_related.generateReport(folderName) # outputs final report with amorce list
		else: 
			print("something went wrong with output files!")
if verbose :
	print("these files have been created :\n",genomicReport,"\n",cdnaReport, "\nhave a nice day !")
sys.exit()





