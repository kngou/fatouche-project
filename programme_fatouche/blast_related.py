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
"""
command_line = input()
/bin/vikings -input eggs.txt -output "spam spam.txt" -cmd "echo '$MONEY'"
>>> args = shlex.split(command_line)
>>> print(args)
['/bin/vikings', '-input', 'eggs.txt', '-output', 'spam spam.txt', '-cmd', "echo '$MONEY'"]
>>> p = subprocess.Popen(args) # Success!

"""


def runBlastN(input_file,subject): 
# query is the longest file it is used as reference 
# subject is the shortest only as convention 
	if not os.path.isfile(subject+".tfa"): # test if the subject file database exists and if it does not it creates it 
		command_line = "makeblastdb -in "+subject+" -dbtype nucl -out " +subject+".tfa"
		cline =  shlex.split(command_line)
		subprocess.Popen(cline)	 
	resultfile = "./result/"+subject+"_result.txt"
	command_line = "blastn -query "+input_file+" -out "+resultfile+" -db "+subject+".tfa"+" -outfmt 7"
	cline =  shlex.split(command_line)
	subprocess.Popen(cline)
	return resultfile

# blastx -query input_file -out output_file
#

# cline = NcbiblastxCommandline._Ncbiblast2SeqCommandline(query="m_cold.fasta", db="nr", evalue=0.001)
# >blastn -query TraesCS5B01G285800.1.fasta -out test2.txt
 # -db cdnaTraesCS5B01G285800.1.tfa
# print (cline)
input_file = "cdnaTemplates.fa"
subject = "genomicTemplates.fa.fasta"

runBlastN(input_file,subject)