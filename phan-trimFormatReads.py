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
import Bio
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from optparse import OptionParser
from Bio import SeqIO
from Bio import SeqRecord
from Bio.Alphabet import IUPAC
#from Bio.Alphabet.IUPAC import ExtendedIUPACProtein


usage= """
This script will take parsed blast data (using parsed info from blast_parse_xml.py) \n will take as input the file "...best.hits" and deduplicate that file\n by taking the read with the lowest e-value score and trimming that read \n according to its proper reading frame and returning the revese complement if necessary
usage: %prog [-g STR] [-f FILE] [-m INT] [-e STR] [-s INT]"""
parser = OptionParser(usage=usage, version="%prog 0.1")

parser.add_option("-g",  dest="gene",
                  help="Specify the gene name which will be used to label output files.")
parser.add_option("-f", dest="blastOutput",
                  help="Specify the input file that is parsed from blast output... *.hits file.")
parser.add_option("-m", dest="minLen",
                  help="Specify the minimum length of reads that will be passed on to rest of analysis, if no min specified all reads will be passed.")
parser.add_option("-e", dest="eValue",
                  help="(OPTIONAL) Evalue cutoff all reads must be =< this evalue to continue,  if no max specified all reads will be considered. Enter in scientific notation (1e-5).")
parser.add_option("-s", dest="score",
                  help="(OPTIONAL) All reads must be >= this score to continue, if no min specified all reads will be considered.")

(options, args) = parser.parse_args()
#Makes sure all mandatory options appear
mandatories = ["gene", "blastOutput"]
for m in mandatories:
	if not options.__dict__[m]:
		print("A mandatory option is missing!\n")
		parser.print_help()
		exit(-1)

#Opens the blast output info file and the new output file
bestParse = os.path.abspath(options.blastOutput)
blast_file = open(bestParse, 'r')

outputFileName = options.gene + "_formatted_query_seqs_NT.fasta" #Can be used for PaPaRa 
outputFile = open(outputFileName, "w")

infoFileName = options.gene + "_query_reads_processing_info.log"
infoFile = open(infoFileName, "w")

ambiguousNucleotideFileName = options.gene + "_unique_reads_with_ambiguous_nucleotides.fasta"
ambiguousFile = open(ambiguousNucleotideFileName, "w")

readsWithStopsFileName = options.gene + "_AA_translations_with_Stops.fasta"
readsWithStopsFile = open(readsWithStopsFileName, "w")

cleanAminoAcidReadsFileName = options.gene + "_AA_translation_stops_removed.fasta" #Can be used for PaPaRa
cleanAminoAcidReadsFile = open(cleanAminoAcidReadsFileName, "w")

dictionaryOutputFileName = options.gene + "_info_DICTIONARY_FILE_for_back_translate.csv"
dictionaryOutputFile = open(dictionaryOutputFileName, "w")

minLength = int(options.minLen)


#Trims lengths according to best blast hit
def endTrim(accession, frame, sequence, queryFrom, hitTo, hitFrom, queryLen, hitLen):

	queryFromNucBase = (queryFrom - 1) * 3
	hitFromPosition = hitFrom - 1

	queryHitPosition = queryFromNucBase + (hitLen - hitFromPosition)
	queryNucLen = queryLen * 3
	
	
	if queryFromNucBase < hitFromPosition: #Trim the front overhang of read
		startTrim = hitFromPosition - queryFromNucBase
		div3 = (hitLen - (frame - 1)) % 3
		endTrim = hitLen - div3
		trimRead = sequence[startTrim:endTrim]


	elif queryHitPosition > queryNucLen: #Trim the end of an overhang read
		maxQueryLen = ((queryLen) - queryFrom) * 3
		endTrim = (hitFrom + maxQueryLen) - 1 

		if frame == 1:
			hitStart = 0
		elif frame == 2:
			hitStart = 1
		elif frame == 3:
			hitStart = 2
		trimRead = sequence[hitStart:endTrim]


	else:
		startTrim = frame - 1
		div3 = (hitLen - (frame - 1)) % 3
		endTrim = hitLen - div3
		trimRead = sequence[startTrim:endTrim]

		

	if len(trimRead) % 3 != 0:
		print("****** NOT div by 3, error in trimming ", accession, " trimRead Len:", len(trimRead), "TrimRead%3:", len(trimRead) % 3)
		print(trimRead)
	return trimRead

#Reads the first few lines of the blast output that contains headers etc...
i = 1
while i < 10:
	blast_file.readline()
	i = i + 1


best_reads = {} #Creates a dictionary, key = accession ID & value = list of read info & read sequence

indexOriginalCount = 0 

for line in blast_file:
		
	if line.startswith(">"):
		indexOriginalCount = indexOriginalCount + 1
		#split the blast info into a line - tab delimited
		line = line.strip("\n")
		info = line.split("\t")

		accession = info[0]
		accession = accession.replace(":", "_").replace("/","_").replace(" ","")
		evalue = float(info[3])
		score = float(info[2])
		frame = int(info[17])
		queryLen = int(info[5])
		queryFrom = int(info[12])
		queryTo = int(info[13])
		hitLen = int(info[7])
		hitFrom = int(info[14])
		hitTo = int(info[15])		

		continue
#Takes care of the breaks in the file between Seed sequences
	elif len(line) < 1:
		continue
	elif line.startswith("#"):
		continue

	elif len(line) > 6:
		read = line.strip("\n")
		
		if frame < 0:
			readToRevComp = Seq(read, generic_dna)
			revComp = str(readToRevComp.reverse_complement())
			sequence = revComp
			frame = abs(frame)
			hitFromOriginal = hitFrom
			hitFrom = (hitLen + 1) - hitTo
			hitTo = hitLen - (hitFromOriginal - 1)

		else:
			sequence = read

		fixedRead = endTrim(accession, frame, sequence, queryFrom, hitTo, hitFrom, queryLen, hitLen)
		

		readInfo = [evalue, frame, read, fixedRead, score]


		test = accession in best_reads #Creates a test to see if the accession is already in the dictionary, will replace if new sequence has better e-value

		if test is False:
			best_reads[accession] = readInfo

		elif test is True:
			old_data = best_reads[accession]
			if evalue < float(old_data[0]):
				del best_reads[accession]
				best_reads[accession] = readInfo

#Reports the orginal number of reads prior to deduplication
infoFile.write("The original number of blast reads prior to deduplication = "+ str(indexOriginalCount) + "\n\n")

#To find original number of unique reads
allKeys = list(best_reads.keys())
numKeys = len(allKeys)
infoFile.write("The original number of unique reads recruited by blast = " + str(numKeys) + "\n\n")

if options.eValue:
	eValue = float(options.eValue)
	allKeys = list(best_reads.keys())
	infoFile.write("e-value cutoff of " + str(eValue) + " used; all reads reported are <= this e-value threshold. \n")
	
	eindex = 0
	for element in allKeys:
		data = best_reads[element]
		if  eValue < float(data[0]):
			del best_reads[element]
			eindex+=1
	infoFile.write("Number of reads which did not make this e-value threshold and were discarded: " + str(eindex) + "\n\n")


if options.score:
	setScore = float(options.score)
	allKeys = list(best_reads.keys())
	infoFile.write("score threshold of " + str(setScore) + " used; all reads reported have score >= this score \n")
	
	sindex = 0
	for element in allKeys:
		data = best_reads[element]
		if setScore > float(data[4]):
			del best_reads[element]
			sindex+=1
	infoFile.write("Number of reads which did not make this score threshold and were discarded: " + str(sindex) + "\n\n")
		

#This sorts all the unique entries in the dictionary and outputs the formmated reads to a file in fasta format
##This also searches for sequences with ambiguous nucleotides and outputs them to a separate file
##This will also add an amino acid translation to the clean unique read dictionary record {best_reads}

ambiguousNucleotides = {}

allAccessions = list(best_reads.keys())
readsWithStops = []


uniqueReadsCount = 0
ambiguousNucleotideCount = 0
cleanUniqueCount = 0
shortReadsCount = 0

for element in allAccessions:
	uniqueReadsCount = uniqueReadsCount + 1
	info = best_reads[element]
	sequence = info[3]
	#removes sequences with ambiguous nucleotides
	if sequence.count("N") > 0 or sequence.count("n") > 0: 
		ambiguousNucleotideCount = ambiguousNucleotideCount + 1
		ambiguousNucleotides[element] = sequence
		del best_reads[element]
	elif len(sequence) <= minLength:
		shortReadsCount = shortReadsCount + 1
		del best_reads[element]
	else:
		cleanUniqueCount = cleanUniqueCount + 1
		outputFile.write(element)
		outputFile.write("\n")
		outputFile.write(sequence)
		outputFile.write("\n")
		coding_DNA = Seq(sequence, IUPAC.unambiguous_dna)
		#coding_DNA = Seq(sequence, ExtendedIUPACProtein)
		AA_seq = str(coding_DNA.translate())
		best_reads[element].append(AA_seq)


		if AA_seq.count("*") > 0:
			readsWithStops.append(element)
			withoutStop = AA_seq.replace("*", "")
			best_reads[element].append(withoutStop)
		else:
			best_reads[element].append(AA_seq)

##prints out a file that has all the reads with stop codons
readsWithStopsCount = 0
for element in readsWithStops:

	readsWithStopsCount = readsWithStopsCount + 1
	readsWithStopsFile.write(element)
	readsWithStopsFile.write("\n")
	##### Location of the read with the stop codons *
	readsWithStopsFile.write(str(best_reads[element][5]))
	readsWithStopsFile.write("\n")

readsWithStopsFile.close()

##prints out a file that has all the unique reads translated into AA space with removed stop codons
goodAccessions = list(best_reads.keys())

dictionaryOutputFile.write("Accession, e-value, abs_frame, raw_read, fixed_read, AA_with_stop, AA_removed_stop \n")

for element in goodAccessions:

	cleanAminoAcidReadsFile.write(element)
	cleanAminoAcidReadsFile.write("\n")
	cleanAminoAcidReadsFile.write(str(best_reads[element][6]))
	cleanAminoAcidReadsFile.write("\n")
	dictionaryOutputFile.write(element)
	values = (best_reads[element])
	for item in values:
		dictionaryOutputFile.write(", ")
		dictionaryOutputFile.write(str(item))

	dictionaryOutputFile.write("\n")


cleanAminoAcidReadsFile.close()	

outputFile.close()

ambiguousAccessions = list(ambiguousNucleotides.keys())

for element in ambiguousAccessions:
	read = ambiguousNucleotides[element]
	ambiguousFile.write(element)
	ambiguousFile.write("\n")
	ambiguousFile.write(read)
	ambiguousFile.write("\n")

ambiguousFile.close()

infoFile.write("The number of unique reads recruited from blast that also met thresholds specified = " + str(uniqueReadsCount) + "\n\n")
infoFile.write("The number of unique reads containing ambiguous nucleotides = " + str(ambiguousNucleotideCount) + "\n")
ambiguousPercent = "%.2f" % ((float(ambiguousNucleotideCount)/uniqueReadsCount)*100)
infoFile.write("The percent of unique reads with ambiguous nucleotides that were removed from analysis = " + str(ambiguousPercent) + "%\n" + "\n")
infoFile.write("The number of unique reads that were less than " + str(minLength) + " bases & removed from analysis = " + str(shortReadsCount) + "\n")
shortReadsPercent = "%.2f" % ((float(shortReadsCount)/uniqueReadsCount)*100)
infoFile.write("The percent of reads that were too short & removed from analysis = " + str(shortReadsPercent) + "%\n" + "\n")
infoFile.write("The number of clean unique reads without ambiguous nucleotides = " + str(cleanUniqueCount) + "\n\n")
infoFile.write("The number of reads with at least one stop codon = " + str(readsWithStopsCount) + "\n")
stopsCalc = "%.2f" % ((float(readsWithStopsCount)/uniqueReadsCount)*100)
infoFile.write("Percent of unique reads with stop codons = " + str(stopsCalc) + "%\n")

blast_file.close()
infoFile.close()