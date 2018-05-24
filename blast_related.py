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
# from Bio.Blast.Applications import NcbiblastxCommandline
import shlex,subprocess


"""
 regroup all functions related to blast 
 and to parse blast result
 it has been decided to output result in text format

"""
def addYourBlastDatabase(subject):
	#create a database file to fit blast command line it return the file in desire format
	subjectDatabaseName =""
	commandLine = "makeblastdb -in "+subject+" -dbtype nucl -out " +subject+".tfa"
	cLine =  shlex.split(commandLine)
	try:		
		process = subprocess.run(commandLine)	
		subjectDatabaseName = ".\\"+subject+".tfa"
		return subjectDatabaseName
	except  :
		print ("Error in database creation. End of the program.")
		exit()


def runBlastN(query,subject): 
	# run blast in command line and generated a result file in csv text format 	
	resultfile = ".\\result\\"+subject+"_result.txt"	 	 	
	subjectfile = ".\\"+subject+".tfa"	 	 	
	#run the blast command line 
	addYourBlastDatabase(subject)		
	command_line = "blastn -query "+query+" -out "+resultfile+" -db "+subjectfile+" -outfmt 7"
	cline =  shlex.split(command_line)
	# print(cline)
	# exit()
	try :
		subprocess.run(command_line)
	except  :	
		subjectfile = addYourBlastDatabase(subject) # try creation of local blast database  but I don't succeed 
		#run the blast command line 		
		command_line = "blastn -query "+query+" -out "+resultfile+" -db "+subjectfile+" -outfmt 7"
		cline =  shlex.split(command_line)
	return resultfile

def parseBlast(resultFile):
	name = ""
	arrayReport = []	
	blastRes = {}	
	try :
		file = open(resultFile,"r")
	except FileNotFoundError:
		print("Your file",resultFile," does not exist !")
		exit()		
	for line in file:
		line =line.strip()	
		if line.startswith('# Query'):
			name = line.replace ('# Query: ','')
			blastRes[name] = arrayReport			
		if not line.startswith("# "):
			arrayReport.append(line)		
		file.close
	return (blastRes)

