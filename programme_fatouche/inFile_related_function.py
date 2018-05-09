#!/usr/bin/env python3
# -*",-coding:Latin-1 -utf8",
import os
import re
import sys
from decimal import *
import FastaSeq


"""
this module groups all functions that deal with fasta inFile

"""



def createResultFolder(folderName="result"):	
	mypath =os.getcwd()+"\\"+folderName # creation of result file in the current directory if it does not exist
	if not os.path.isdir(mypath):
		os.makedirs(mypath)
		print( "A new folder has been created.",mypath)
	return mypath


def ExtractFasta(filename):

	"""
		extract sequence from infile in fasta format and put it on a dictionnary
		return a dictionnary
	"""
	FastaSeq={}
	path = os.getcwd()
	fileIN = path+"\\"+filename	
	try :
		file = open(fileIN,"r")
	except FileNotFoundError:
		print("Your last input file does not exist ! Check your file name")
		exit()	
	seq =""	
	name =""
	for line in file:
		line =line.strip()			
		if line.startswith('>'):
			name = line.replace(">","")
			seq=""
		else:
			seq += line
			FastaSeq[name] = seq
	file.close	
	return FastaSeq

def SequenceFromExonFile(filename):
	
	ExonSeq={}	
	fileIN = ".\\result\\"+filename+"_"+"exon_list.txt"
	# print (fileIN)
	try :
		file = open(fileIN,"r")
	except FileNotFoundError:
		print("Your last input file does not exist ! Check your file name")
		exit()	
	seq =""	
	name =""
	for line in file:
		line =line.strip()			
		if line.startswith('exon'):
			name = line
			seq=""
		else:
			seq += line
			ExonSeq[name] = seq
	file.close	
	return ExonSeq

def isExist(exonSelected,filename):
	ExonSeq = SequenceFromExonFile(filename)
	if ExonSeq.get(exonSelected):
		return True
	else:
		return False

def getExon(exonSelected,filename):
	ExonSeq = SequenceFromExonFile(filename)
	return ExonSeq.get(exonSelected)
	



