#!/usr/bin/env python

#This file is part of the PHANS-C pipeline.

#the PHANS-C pipeline is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#The PHANS-C pipeline is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with the PHANS-C pipeline.  If not, see <http://www.gnu.org/licenses/>.

#Copyright: University of Washington 2017
#Author: Jaclyn K Saunders

import sys
import os
import re
import os.path
import shutil
from optparse import OptionParser

usage= """
Takes an alignment file from PaPaRa -- identifies paired end reads and merges these reads for placement with EPA

usage: %prog [-p FILE] [-r INT] [-o FILE] [-l INT] """

parser = OptionParser(usage=usage, version="%prog 0.1")

parser.add_option("-p", "--PaPaRa file", dest="alnFile",
                  help="Specify the output file from PaPaRa that has the reference aln and metagenome reads aln combined",
                  metavar="FILE")
parser.add_option("-r", "--number of reference sequences in integer format", dest="refInt",
                  help="Specify the number of reference sequences in the alignment (# of sequences in reference tree)",
                  metavar="INT")
parser.add_option("-o", "--output alingment filename", dest="outputFile",
                  help="Specify the output file name",
                  metavar="FILE")
parser.add_option("-l", "--minimum length of sequence", dest="minLen",
				  help="Specify the minimum number of residues that must be in each read for placement; default set to 45 bases/residues",
				  metavar="INT")

(options, args) = parser.parse_args()

#Makes sure all mandatory options appear
mandatories = ["alnFile", "refInt"]
for m in mandatories:
	if not options.__dict__[m]:
		print("A mandatory option is missing!\n See the HELP menu - 'spaceJoin.py -h'") 
		parser.print_help()
		exit(-1)
#Sets the default outputfile name	
if options.outputFile == None:
	outputName = options.alnFile
	outputFileName = outputName + "_joinedReadsForPlacement.phymlAln"
else: 
	outputName = options.outputFile
	outputFileName = outputName + "_joinedReadsForPlacement.phymlAln"

#Sets default minimum length to 45 if none specified
if options.minLen == None:
	minLength = 45
else: 
	minLength = int(options.minLen)

#Set variables based on collected input
paparaAlnFileName = os.path.abspath(options.alnFile)
referenceSeqCount = options.refInt
refCount = int(referenceSeqCount)

outputFile = open(outputFileName, "w")
(outputFileNamePre, outputFileNameSuf) = os.path.splitext(outputFileName)


logFileName = "Log_spaceJoin_" + outputFileNamePre + ".txt"
logFile = open(logFileName, "w")

alnFileOpen = open(paparaAlnFileName, 'r')

alnFileFirstLine = alnFileOpen.readline()

splitFirstLine = alnFileFirstLine.split(" ")

numSeqsInFile = splitFirstLine[0]
numPositionsInAln = splitFirstLine[1]

###To test is a string represents and int value
def RepresentsInt(string):
	try:
		int(string)
		return True
	except ValueError:
		return False

###Generate list to output later
ouputLines = []

###Create dictionary of reads for merging later
reads = {}
preJoinedReads = {}

###Read and open the rest of the aln file
referenceSequences = []

index = 1
for line in alnFileOpen:
	
	if index <= refCount: #Removes the reference sequences
		referenceSequences.append(line)
		outputFile.write(line)
		index += 1
		
	else: #Metagenomic reads in alignment
		numSpaces = line.count(" ")
		line = line.replace("\n","")
		lineSplit = line.split(" ")
		readID = lineSplit[0]
		readSequence = lineSplit[numSpaces]
		readKey = readID[0:-1]
		readPair = readID[-1]
		value = []
#To handle the new combined reads in the blastdb which are pre-combined and lack an _1 or _2; may have _0 to signify pre-joined
		if readID[-2:-1] != "_": #checks to make sure a true read pair at the end and not just an int at the end of the read name
			value = [numSpaces, readID, readSequence, "PreJoined"]
		else:
			if RepresentsInt(readPair) is False:
				value = [numSpaces, readID, readSequence, "PreJoined"]

			else:
				readPair = int(readPair)
				if readPair == 1:
					value = [numSpaces, readID, readSequence, 0, 0]
				elif readPair == 2:
					value = [numSpaces, 0, 0, readID, readSequence]
				elif readPair == 0:
					value = [numSpaces, readID, readSequence, "PreJoined"]
				else:
					print("Error #1")
					print("\t" + readID) 

		test = readKey in reads

		if test is False:
			### Add data
			reads[readKey] = value

		elif test is True:
			### Pull old data & update
			if readPair == 1:
				oldValue = reads[readKey]
				read2ID = oldValue[3]
				read2Seq = oldValue[4]
				newValue = [numSpaces, readID, readSequence, read2ID, read2Seq]
				del reads[readKey]
				reads[readKey] = newValue
			elif readPair == 2:
				oldValue = reads[readKey]
				read1ID = oldValue[1]
				read1Seq = oldValue[2]
				newValue = [numSpaces, read1ID, read1Seq, readID, readSequence]
				del reads[readKey]
				reads[readKey] = newValue
			else:
				print("Error #2")
				print("\t\t" + readID + "\t" + readPair)

allKeys = list(reads.keys())

preJoinedReads = [] ###To handle the sequences without _1 or _2; may have _0 to signify pre-joined
bothReads = []
onlyRead1 = []
onlyRead2 = []

for element in allKeys:
	values = reads[element]
	if values[3] == "PreJoined":
		preJoinedReads.append(element)
	else:

		if type(values[1]) == str:
			read1 = True
		else:
			read1 = False

		if type(values[3]) == str:
			read2 = True
		else:
			read2 = False

		if read1 == True and read2 == True:
			bothReads.append(element)
		elif read1 == True and read2 == False:
			onlyRead1.append(element)
		elif read1 == False and read2 == True:
			onlyRead2.append(element)
		else:
			print("Error #3 \t \t" + element)

numPreJoinedReads = len(preJoinedReads)
numBothReads = len(bothReads)
numOnlyRead1 = len(onlyRead1)
numOnlyRead2 = len(onlyRead2)


###Whole dataset
logFile.write("The following of the input alignment are the following: \n\n")
logFile.write("Pre-Joind BOTH READS: " + str(numPreJoinedReads) + "\n") 
logFile.write("BOTH READS: " + str(numBothReads) + "\n") 
logFile.write("ONLY READ 1: " + str(numOnlyRead1) + "\n") 
logFile.write("ONLY READ 2: " + str(numOnlyRead2) + "\n") 
logFile.write("\n\n")

#######Create a definition for carrying out length criteria#########
def lengthCriteria(sequence):
	sequenceNoGaps = sequence.replace("-","")
	sequenceLen = len(sequenceNoGaps)
	#####Can alter length criteria here######
	if sequenceLen >= minLength:
		return True
	else:
		return False

longPreJoinedReads = 0
longRead1 = 0
longRead2 = 0
longBothReads = 0

for element in preJoinedReads:
	data = reads[element]
	seqName = data[1]
	spaces= data[0]
	alnSeq = data[2]
	numSpaces = " "*spaces
	newLine = "\n"
	outputToFile = seqName + numSpaces + alnSeq + newLine
	if lengthCriteria(alnSeq) == True:
		outputFile.write(outputToFile)
		longPreJoinedReads+=1

for element in onlyRead1:
	data = reads[element]
	seqName = data[1]
	spaces = data[0]
	alnSeq = data[2]
	numSpaces = " "*spaces
	newLine = "\n"
	outputToFile = seqName + numSpaces + alnSeq + newLine
	if lengthCriteria(alnSeq) == True:
		outputFile.write(outputToFile)
		longRead1+=1
	#else:
		#print "SEQUENCE " + seqName + " is too short!"

for element in onlyRead2:
	data = reads[element]
	seqName = data[3]
	spaces = data[0]
	alnSeq = data[4]
	numSpaces = " "*spaces
	newLine = "\n" 
	outputToFile = seqName + numSpaces + alnSeq + newLine
	if lengthCriteria(alnSeq) == True:
		outputFile.write(outputToFile)
		longRead2+=1
	#else:
		#print "SEQUENCE " + seqName + " is too short!"

for element in bothReads:
	data = reads[element]
	seqName = element + "3"
	spaces = data[0]
	numSpaces = " "*spaces
	read1 = data[2]
	read2 = data[4]
	lenSeq = len(read1)
	index = 0
	joinedRead = []

	while lenSeq > index:
		positionRead1 = read1[index]
		positionRead2 = read2[index]
		if positionRead1 == positionRead2:
			joinedRead.append(positionRead1)
		elif positionRead1 != positionRead2:
			if positionRead1 != "-":
				joinedRead.append(positionRead1)
			else:
				joinedRead.append(positionRead2)
		index += 1

	stringRead = ''.join(joinedRead)
	outputToFile = seqName + numSpaces + stringRead + newLine

	if lengthCriteria(stringRead) == True:
		outputFile.write(outputToFile)
		longBothReads+=1


logFile.write("*****************************************************************************************\n")
logFile.write("The following are the set of reads in the output file that are longer than " + str(minLength) + " \n")
logFile.write("LONG Pre-Joined BOTH READS: " + str(longPreJoinedReads) + "\n") 
logFile.write("LONG BOTH READS: " + str(longBothReads) + "\n") 
logFile.write("LONG ONLY READ 1: " + str(longRead1) + "\n")
logFile.write("LONG ONLY READ 2: " + str(longRead2) + "\n") 
numAllReads = numBothReads*2 + numOnlyRead1 + numOnlyRead2 + numPreJoinedReads*2
numLongReads = longBothReads*2 + longRead1 + longRead2 + longPreJoinedReads*2
percentTotal = '%.2f' % ((float(numLongReads) / numAllReads)*100)
logFile.write("\nPerecent of reads from total dataset that are longer than " + str(minLength) + ": " + str(percentTotal) + "%\n")
percentPairedEndReads = '%.2f' % (((float(longBothReads*2) + float(longPreJoinedReads*2)) / (longBothReads*2 + longRead1 + longRead2 + longPreJoinedReads*2))*100)
percentPairedEndReads = str(percentPairedEndReads)
logFile.write("\nPercent of long reads with pairs: " + percentPairedEndReads + "%\n")
logFile.write("\n\n")
logFile.write("*****************************************************************************************\n")
logFile.write("Reads that are written to output file \n")
logFile.write("\nReads with mate pairs joined: \n")
logFile.write("3, ".join(bothReads))
logFile.write("\n\nReads with only Read 1: \n")
logFile.write("1, ".join(onlyRead1))
logFile.write("\n\nReads with only Read 2: \n")
logFile.write("2, ".join(onlyRead2))
logFile.write("\n\nReads that were already pre-Joined: \n")
logFile.write("0, ".join(preJoinedReads))

###To write header in new alignment file:
numberReadsToPlace = longBothReads + longRead1 + longRead2 + longPreJoinedReads
totalSequencesInFile = numberReadsToPlace + refCount

outputFile.close()

###Depending on size of file and amount of memory available, this may not work -- in which case you would need to read line by line
tempFile = open(outputFileName, "r")
tempRead = tempFile.read()
tempFile.close()

finalFile = open(outputFileName, "w")
finalFile.write(str(totalSequencesInFile) + " " + numPositionsInAln)
finalFile.write(tempRead)
finalFile.close()

outputFile.close()
logFile.close()
alnFileOpen.close()