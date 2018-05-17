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
from Bio.Blast.Applications import NcbiblastxCommandline
import shlex,subprocess


"""
 regroup all functions related to blast command line
 and to parse blast result
 it has been decided to output result in text format

"""
def addYourBlastDatabase(subject):
	#create a database file to fit blast command line it return the file in desire format
	subjectDatabaseName =""
	commandLine = "makeblastdb -in "+subject+" -dbtype nucl -out " +subject+".tfa"
	cLine =  shlex.split(commandLine)
	subprocess.Popen(cLine)	
	subjectDatabaseName = subject+".tfa"
	return subjectDatabaseName



def runBlastN(input_file,subject): 
	# run blast in command line and generated a result file in csv text format 	
	resultfile = "./result/"+subject+"_result.txt"	 
	 	# create the subject file database 
	subjectfile = subject+".tfa"
	commandLine = "makeblastdb -in "+subject+" -dbtype nucl -out " +subject+".tfa"
	cLine =  shlex.split(commandLine)
	subprocess.Popen(cLine)
		#run the blast command line 		
	command_line = "blastn -query "+input_file+" -out "+resultfile+" -db "+subjectfile+" -outfmt 7"
	cline =  shlex.split(command_line)
	subprocess.Popen(cline)
	return resultfile

def parseBlast(resultFile):
	try :
		file = open(resultFile,"r")
	except FileNotFoundError:
		print("Your file",resultFile," does not exist !")
		exit()	
	
	for line in file:
		line =line.strip()	
		array = line.split(",")
		# print (array)
		if not line.startswith('#'):
			queryAccVer, subjectAccVer,identity, alignmentLength, mismatches, gapOpens, qStart, qEnd, sStart, sEnd, evalue, bitScore = line.split('\t')
		file.close	
		# print (queryAccVer)


input_file = "cdnaTemplates.fa"
subject = "genomicTemplates.fa"

results = runBlastN(input_file,subject)
parseBlast(results)






