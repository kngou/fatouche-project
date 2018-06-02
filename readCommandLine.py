#!/usr/bin/python3
# -*",-coding:Latin-1 -utf8",
import os
import re
import sys
import argparse 



def commandLine():
	folderName = ""
	master = ""
	genomicFile = ""
	cdnaFile = ""
	verbose = ""
	
	parser = argparse.ArgumentParser()

	parser.add_argument("-m", "--master", action ="store_true",default = False, dest= "master", 
						help = "It searches on the current folder all files with fasta extension and it uses it as values for genomic file name. Use it if your current folder contains only the fasta files that you need. Pay attention, that the corresponding cdna file must have fa as extension unless you will retrieve an error message")
	parser.add_argument("-g", "--genomicFile ",
	                   action ="store", type =str, dest= "genomicFile", nargs ="*",	                   
	                   help ="genomic file name to input. If you input a list of files, file names must be separated by a coma.")
	parser.add_argument("-v", "--verbose",
	                   action ="store_true", dest="verbose", default = False,
	                   help ="choose this option if you want to output a message for each step.")
	parser.add_argument("-c", "--cdnaFile ",
	                  action ="store", type =str, dest = "cdnaFile", nargs ="*",
	                  help ="cdna file name list to input. Each file name must be separated by a coma. Files must be in the same order as its corresponding genomic file, or it must have the same name.")
	parser.add_argument("-f", "--folderName ",action ="store", type =str, dest = "folderName", default="result",
					  help ="folder name for ouput files")
	options = parser.parse_args()
	master = options.master
	verbose = options.verbose
	folderName = options.folderName
	# -- input file names are stored on list each file is separated by a coma 
	genomicList = options.genomicFile
	cdnaList = options.cdnaFile	
	return folderName,master,genomicList,cdnaList,verbose
