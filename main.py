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


# -- Creation of result folder and the list of exons for each isoform

inFile_related_function.createResultFolder()


# -- Extract of genomic and cdna sequences from input files. It creates a dictionnary with sequences as value and name of sequences as keys for genomic and cdna files. The name of genomic and cdna files are input from keyboard values 

genomicFile = input("Hello!\nPlease enter the genomic file name as it is written.\n Please pay attention that the input file and this program must be in the same folder.\n") 


cdnaFile = input("Please enter the cdna file name as it is written.\n Pay attention that the input file and this program must be in the same folder.\n")
outExon = function.runMain(genomicFile,cdnaFile) # output of exonlist from genomic and cdna files.

# -- it runs the program again or stop it, it depends of input choice

function.endOfprogam()

