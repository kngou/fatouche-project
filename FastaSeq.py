#!/usr/bin/python3.4
# -*",-coding:Latin-1 -utf8",
import os
import re
import sys
import function

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
	GCrateFW = function.GCcontent(forward)
	TmrateFW = function.TmCheck(forward)
	TmrateRV = function.TmCheck(reverse)
	GCrateRV = function.GCcontent(reverse)
	if GCrateFW <= 40 or GCrateFW >= 60:
		forward = "not available '%'GC issues",GCrateFW
	if GCrateRV <= 40 or GCrateRV >= 60:
		reverse = "not available '%'GC issues",GCrateRV
	if TmrateFW <= 50 or TmrateFW >= 61:
		forward = "not available Tm issues",TmrateFW
	if TmrateRV <= 50 or TmrateRV >= 61:
		reverse = "not available Tm issues",TmrateRV	
	return (forward,reverse)
	if __name__ == "__main__":
	    import doctest
	    doctest.testmod()


def genomicPrimers(sequence,cuttOff=20) :
	forward = "not available "
	reverse = "not available "
	if len(sequence)>= 40:
		forward = sequence[0:cuttOff]
		reverseSeq = RevCom(sequence)
		reverse = reverseSeq[10:cuttOff+10]		
		GCrateFW = function.GCcontent(forward)
		TmrateFW = function.TmCheck(forward)
		TmrateRV = function.TmCheck(reverse)
		GCrateRV = function.GCcontent(reverse)
		if GCrateFW <= 40 or GCrateFW >= 60:
			forward = "not available '%'GC issues",GCrateFW
		if GCrateRV <= 40 or GCrateRV >= 60:
			reverse = "not available '%'GC issues",GCrateRV
		if TmrateFW <= 50 or TmrateFW >= 61:
			forward = "not available Tm issues",TmrateFW
		if TmrateRV <= 50 or TmrateRV >= 61:
			reverse = "not available Tm issues",TmrateRV	
	else: 
		print("check your sequence !")
		
	return (forward,reverse)
	if __name__ == "__main__":
	    import doctest
	    doctest.testmod()