#!/usr/bin/python3.4
# -*",-coding:Latin-1 -utf8",
import os
import re
import sys


"""
This module is used to group function dealing with sequence manipulation
it takes the input sequence and return the reverse and the complement of the sequence
input sequence must be a string in uppercase.

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


def getPrimers(sequence,cuttOff=20):
	"""
	getPrimers function 
	-------------------
	From input sequence it outputs forwards and reverse amorces. it has been decided that the length of the amorces is 20 nucleotides, but it can be changed by adding an integer as second option, when this function is called.
	Pay attention that the sequence is on 5'- 3' so as the two amorces. If you use serial cloner application you do not have to reverse the second amorce.
	:param sequence: The sequence that you want to extract the two amorces.
	:param cuttOff: define the length of the amorce. It is an optional parameter, it is set on 20 by defaut.
	:type sequence: string
	:type cuttOff: int 
	return: Two strings the reverse and forward amorces

	:Example:

		sequence = "ATGAAGATGAAGATTATGAAGATGAAGATTAAGATGAAGATGAAGATTATGAAGATGAAGATTAAGAAGATGAAGATGAAGATTATGAAGATGAAGATTAAGATGAAGATGAAGATTATGAAGATGAAGATTAAGAAGCTGA"
	 >>> getPrimers(sequence,cuttOff=20)
		forward = "ATGAAGATGAAGATTATGAA"
		reverse = "TCAGCTTCTTAATCTTCATC"

	""" 
	forward = sequence[0:cuttOff]
	reverseSeq = RevCom(sequence)
	reverse = reverseSeq[0:cuttOff]	
	return (forward,reverse)
	if __name__ == "__main__":
	    import doctest
	    doctest.testmod()