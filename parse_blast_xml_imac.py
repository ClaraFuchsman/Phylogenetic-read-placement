#!/usr/bin/env python

import sys
import os.path
import subprocess
from optparse import OptionParser
import Bio
from Bio import SeqIO
from Bio.Blast import NCBIXML
import itertools
from Bio.Blast.Record import Alignment


#### Collect Input ####
#######################
usage="""
Takes a file with >=1 XML format blast results and makes a tab delimeted output file.
One hit per line. If no outfile is specified then the outfile is the name of the infile, with "best"
as the suffix.
usage: %prog -f FILE"""

parser = OptionParser(usage=usage, version="%prog 1.2")


parser.add_option("-f", "--in_file", metavar="FILE", dest="in_file",
				  help="Specify the blast XML format FILE that you wish to use as input")

parser.add_option("-o", "--out_file", metavar="FILE", dest="out_file", default = False,
				  help='Specify the outfile that you wish to use. If nothing is specified, the '
					   'name and location of the infile will be used, with a "best" as the suffix.')

parser.add_option("-n", "--number_hits", metavar="INT", default = "1", dest="blast_rows",
				  help="Specify the number of hits you want from each blast record. A 0 "
					   "means return all results. Default is 1.")

parser.add_option("-s", "--sort", default = False, dest="sort", action="store_true",
				  help="Sort by evalue.")

parser.add_option("-l", "--long", default = False, dest="long_results", action="store_true",
				  help="Return extra information from blast results")

parser.add_option("-e", "--extract_hits", default = False, dest="extract_hits",
				  help="Specify db to extract hits from.")

parser.add_option("-c", "--cores", default = '2', dest="num_cores",
				  help="Specify number of cpu cores to use. Default=2")

				  
(options, args) = parser.parse_args()



if not options.in_file:
	print "You must specify an xml input file"
	sys.exit()
	
in_file = os.path.abspath(options.in_file)

(in_file_path, in_file_whole_name) = os.path.split(in_file)
(in_file_base, in_file_ext) = os.path.splitext(in_file_whole_name)

if options.out_file:
	out_file = os.path.abspath(options.out_file)
else:
	out_file = os.path.join(in_file_path, in_file_base + ".best")


extract_hits = options.extract_hits
if extract_hits:
	extract_hits_file = out_file + ".hits"

blast_rows = int(options.blast_rows)

sort = options.sort
long_results = options.long_results
num_cores = int(options.num_cores)


## Functions ##
###############

# This function returns the sortable part from each line. Then we use it
# as part of a key function in the sort.
def by_evalue(line):
	value = line[4] 
	try:
		value = float(value)
	except:
		pass
	return value

def by_query(line):
	value = line[0] 
	try:
		value = float(value)
	except:
		pass
	return value
	

## Classes ##
#############
class Alignment(Bio.Blast.Record.Alignment):
	"""Extends the Bio.Blast.Record.Alignment class to search Genbank to find the likely GI number of a hit.
	Simply returns the first GI if there is more than one. Adds a method that downloads and stores a
	SeqRecord object as an attribute. Also adds a couple of methods to return cleaned titles."""
		
	def __init__(self):
		if not hasattr(self, "proteinGI"):
			self.proteinGI = ""

		if not hasattr(self, "nucleotideGI"):
			self.nucleotideGI = ""

		if not hasattr(self, "proteinFeature"):
			self.proteinFeature = ""

		if not hasattr(self, "nucleotideFeature"):
			self.nucleotideFeature = ""


		if not hasattr(self, "extraCleanTitle"):
			self.extraCleanTitle = self.extraCleanUpTitle()

		if not hasattr(self, "cleanTitle"):
			self.cleanTitle = self.cleanUpTitle()
	
	def convertInstance(self):
		"""Call this function just after converting an instance to this class to do __init__ like things."""
		if not hasattr(self, "proteinGI"):
			self.proteinGI = ""
			
		if not hasattr(self, "nucleotideGI"):
			self.nucleotideGI = ""

		if not hasattr(self, "cleanTitle"):
			self.cleanTitle = self.cleanUpTitle()

		if not hasattr(self, "extraCleanTitle"):
			self.extraCleanTitle = self.extraCleanUpTitle()

		if not hasattr(self, "proteinFeature"):
			self.proteinFeature = ""
		
		if not hasattr(self, "nucleotideFeature"):
			self.nucleotideFeature = ""
			
	
	def getProteinGI(self):
		"""Makes a best guess as to GI number. Stores the number in the "proteinGI" attribute,
		so we don't have to look it up multiple times."""
		
		if not self.proteinGI:
			try:
				giList = GenBank.search_for(self.accession, "protein")
				if len(giList) < 1:
					self.proteinGI = "no_GI_found"	
				else:
					if len(giList) > 1:
						print "found multiple possible GI numbers with %s as a query: %s" % (self.accession, giList)
					self.proteinGI = giList[0]
			except:
				print "Woops! Couldn't get protein feature from NCBI. Maybe network is slow? I'll wait 30 seconds, then try again."
				sleep(30)
				giList = GenBank.search_for(self.accession, database="protein")
				if len(giList) < 1:
					self.proteinGI = "no_GI_found"	
				else:
					if len(giList) > 1:
						print "Found multiple possible GI numbers with %s as a query: %s" % (self.accession, giList)
					self.proteinGI = giList[0]
		return self.proteinGI



	def getNucleotideGI(self):
		"""Makes a best guess as to GI number. Stores the number in the "nucleotideGI" attribute,
		so we don't have to look it up multiple times."""
		
		if not self.nucleotideGI:
			try:
				giList = GenBank.search_for(self.accession, "nucleotide")
				if len(giList) < 1:
					self.nucleotideGI = "no_GI_found"	
				else:
					if len(giList) > 1:
						print "found multiple possible GI numbers with %s as a query: %s" % (self.accession, giList)
					self.nucleotideGI = giList[0]
			except:
				print "Woops! Couldn't get protein feature from NCBI. Maybe network is slow? I'll wait 30 seconds, then try again."
				sleep(30)
				giList = GenBank.search_for(self.accession, database="nucleotide")
				if len(giList) < 1:
					self.nucleotideGI = "no_GI_found"	
				else:
					if len(giList) > 1:
						print "found multiple possible GI numbers with %s as a query: %s" % (self.accession, giList)
					self.nucleotideGI = giList[0]
		return self.nucleotideGI

	
	
	def getProteinFeature(self):
		"""Tries to get whole aa feature from Genbank, but checks first to see if the instance already has it."""
		
		if not self.proteinFeature:
			try:
				record_parser = GenBank.FeatureParser()
				ncbiDict = GenBank.NCBIDictionary('protein', 'genbank', parser = record_parser)
				self.proteinFeature = ncbiDict[self.getProteinGI()]
				self.proteinFeature.__class__ = SeqRecord
				self.proteinFeature.convertInstance()
			except:
				print "Woops! Couldn't get protein feature from NCBI. Maybe network is slow? I'll wait 30 seconds, then try again."
				sleep(30)
				record_parser = GenBank.FeatureParser()
				ncbiDict = GenBank.NCBIDictionary('protein', 'genbank', parser = record_parser)
				self.proteinFeature = ncbiDict[self.getProteinGI()]
				self.proteinFeature.__class__ = SeqRecord
				self.proteinFeature.convertInstance()
		return self.proteinFeature


	def getNucleotideFeature(self):
		"""Tries to get whole nucleotide feature from Genbank, but checks first to see if the instance already has it."""
		
		if not self.nucleotideFeature:
			try:
				record_parser = GenBank.FeatureParser()
				ncbiDict = GenBank.NCBIDictionary('nucleotide', 'genbank', parser = record_parser)
				self.nucleotideFeature = ncbiDict[self.getNucleotideGI()]
				self.nucleotideFeature.__class__ = SeqRecord
				self.nucleotideFeature.convertInstance()
			except:
				print "Woops! Couldn't get nucleotide feature from NCBI. Maybe network is slow? I'll wait 30 seconds, then try again."
				sleep(30)
				record_parser = GenBank.FeatureParser()
				ncbiDict = GenBank.NCBIDictionary('nucleotide', 'genbank', parser = record_parser)
				self.nucleotideFeature = ncbiDict[self.getNucleotideGI()]
				self.nucleotideFeature.__class__ = SeqRecord
				self.nucleotideFeature.convertInstance()

		return self.nucleotideFeature
		
	
	def cleanUpTitle(self):
		"""Make a cleaner version of self.title that is a bit more human readable."""
		clean = self.title.replace("lcl|","")
		clean = clean.split('>', 1)
		clean = clean[0]
		return str(clean)
		
	def extraCleanUpTitle(self):
		"""Make a cleaner version of self.title that is even more human readable. Also splits off position and color from
		our fasta title format."""
		clean = self.cleanTitle
		if clean.count("|") == 4 and clean.startswith("gi") == False:
			clean = clean.split("|")
			clean = clean[0] + "|" + clean[1] + "|" + clean[2]
		return str(clean)




## Program ##
#############


in_file_handle = open(in_file, 'rU')
results_list = []


# Workaround crappy blast behavior. Need to test whether Blast has output plaintext error 
# messages in front of valid XML block If so, need to skip until we get to XML part. Once
# there we set the cursor to the start of the xml block, and proceed as usual
xml_start_line = 0
file_seek_position = 0
max_lines = 1000 #This is the number of lines of error lines we will tolerate before we give up
while xml_start_line < max_lines:
	line = in_file_handle.readline()
	if not line.startswith('<?xml'):
		file_seek_position = in_file_handle.tell()
		xml_start_line = xml_start_line + 1
	else:
		in_file_handle.seek(file_seek_position, 0)
		xml_start_line = max_lines

parsed_blasts = NCBIXML.parse(in_file_handle)


for parsed_blast_record in parsed_blasts:
	if len(parsed_blast_record.alignments) > 0: #Check to see if we have any hits at all.
		i = 0
		while (i < blast_rows or blast_rows == 0) and len(parsed_blast_record.alignments) > i:
			hit = parsed_blast_record.alignments[i]
			hit.__class__ = Alignment
			hit.convertInstance()
			row = "%s	%s	%s	%s	%s	%s\n" % (
					parsed_blast_record.query, hit.accession, hit.cleanTitle, hit.hsps[0].score,
					hit.hsps[0].expect, str(i+1))
			if long_results:
				if hit.hsps[0].gaps == (None, None):
					gaps = 0
				else:
					gaps = hit.hsps[0].gaps
				if len(hit.hsps[0].frame) > 1:
					query_frame = hit.hsps[0].frame[0]
					hit_frame = hit.hsps[0].frame[1]
				else:
					query_frame = 'NA'
					hit_frame = hit_frame = hit.hsps[0].frame[0]
					
				row = "%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s	%s\n" % (
						parsed_blast_record.query, hit.accession, hit.cleanTitle, hit.hsps[0].score,
						hit.hsps[0].expect, str(i+1), parsed_blast_record.query_letters, hit.hit_id,
						hit.length, hit.hsps[0].align_length, hit.hsps[0].identities,
						hit.hsps[0].positives, gaps, hit.hsps[0].query_start, hit.hsps[0].query_end,
						hit.hsps[0].sbjct_start, hit.hsps[0].sbjct_end, query_frame, hit_frame )
			else:
				row = "%s	%s	%s	%s	%s	%s\n" % (
						parsed_blast_record.query, hit.accession, hit.cleanTitle, hit.hsps[0].score,
						hit.hsps[0].expect, str(i+1))
			results_list.append(row)
			i = i + 1
	else: # This is what to do if we didn't have any hits.
		if long_results:	
			row = ("%s	None	None	None	None	None	None	None	None	None	"
				   "None	None	None	None	None	None	None	None	None\n") % (parsed_blast_record.query)
		else:
			row = "%s	None	None	None	None	None\n" % (parsed_blast_record.query)
		results_list.append(row)


#Convert tab delimited string to list, split on tabs. Why did I not make it like that in first place?
temp_list = []
for line in results_list:
	temp_list.append(line.split("	"))
results_list = temp_list
	
if sort:
	results_list.sort(key=by_evalue)


if long_results:
	if extract_hits:
		extract_title_line = ("Each hit has a fasta title in the following tab delimited format:\n"
			"Query	Hit Accession	Description Score	Evalue	Hit Number	"
			"Query Length	Hit ID	Hit Length	Align Length	Identities	Positives	"
			"Gaps	Query From	Query To	Hit From	Hit To	Frame Query Frame Hit\n\n")
		
	title_Line = ("Query	Hit Accession	Description Score	Evalue	Hit Number	"
		"Query Length	Hit ID	Hit Length	Align Length	Identities	Positives	"
		"Gaps	Query From	Query To	Hit From	Hit To	Frame Query Frame Hit\n")
else:
	if extract_hits:
		extract_title_line = ("Each hit has a fasta title in the following tab delimited format:\n"
			"Query	Hit Accession	Description Score	Evalue	Hit Number\n\n")

	title_Line = "Query Hit Accession	Description Score	Evalue	Hit Number\n"
	


out_file_handle = open(out_file, 'w')
out_file_handle.write(title_Line)

if extract_hits:
	extract_hits_tmp_file_handles = []
	i = 0
	while i < num_cores:
		extract_hits_tmp_file_handles.append(open('{}_{}.tmp'.format(extract_hits_file, i), 'w+'))
		i += 1
	extract_hits_tmp_file_handle_cycle = itertools.cycle(extract_hits_tmp_file_handles)



#Not sure if already in sort order or not. 
results_list.sort(key=by_query)


count_kept = 0
count_total = 0
seen_results = set()
for line in results_list:
	count_total += 1
	out_file_handle.write("\t".join(line))
	
	if extract_hits:
		#Distribute all unique hits to tmp files for later use.
		if line[1] not in seen_results and line[1] != 'None':
			if line[1] == 'None':
				print line
			f = extract_hits_tmp_file_handle_cycle.next()
			f.write('{}\n'.format(line[1]))
			seen_results.add(line[1])
			count_kept += 1
del seen_results #This set can be big, so get rid of it immediately.

if extract_hits:
	print "{} of {} were unique hits ({:.2%}). We will extract these from blast db".format(
		count_kept, count_total, float(count_kept)/count_total)			

if extract_hits:
	print 'Extracting fasta hits from the reference blast db. This can take a while.'
	# Start num_cores blastdbcmd processes.
	processes = []
	for f in extract_hits_tmp_file_handles:
		f.close() #close bc an external program will need it
		extract_hits_tmp_file = f.name
		cmd ='/Applications/blast/bin/blastdbcmd -db {} -entry_batch {} -out {}'.format(
			  extract_hits, extract_hits_tmp_file, extract_hits_tmp_file + '.fasta')
		
		processes.append(subprocess.Popen(cmd, shell=True))
	
	#Now wait for them all to finish
	for p in processes:
		status = p.wait()
		if status:
			print 'Oh crap, blastdbcmd failed, returning: {} {}'.format(status, p.stderr)
			sys.exit()
	print 'Finished extracting fasta hits from the reference blast db.'
		

	#Open fastas generated by blastdbcmd and cat them together, then index them.
	with open(extract_hits_file + '_all.fasta', 'w') as ref_handle:
		for f in extract_hits_tmp_file_handles:
			with open(f.name + '.fasta') as tmp_fasta_handle:
				ref_handle.write(tmp_fasta_handle.read())
	fasta_hits = SeqIO.index(ref_handle.name, 'fasta') #ref_handle was implicitly closed

	
	with open(extract_hits_file, 'w') as f:
		f.write(extract_title_line)
		last_query_seen = ''
		for line in results_list:
			query = line[0]
			hit_name = line[1]
			if query != last_query_seen:
				last_query_seen = query
				f.write('\n\n#####\n#Hits below from query {}\n#####\n\n'.format(query))

			if hit_name != 'None':
				try: #Try because stupid blast often auto-prepends 'lcl|' to identifiers
					hit_fasta_string = fasta_hits['lcl|' + hit_name].format('fasta')
				except KeyError:
					hit_fasta_string = fasta_hits[hit_name].format('fasta') 
		
				if not hit_fasta_string[0] == '>':
					print "Uh oh. Something went wrong, we should have gotten a fasta record back: {}".format(hit_fasta_string)
					sys.exit()
		
				full_fasta_name = "\t".join(line[1:])
				fasta_seq = ''.join(hit_fasta_string.split()[1:])
				fasta_record = ">{0}{1}\n".format(full_fasta_name, fasta_seq) #Second item is fasta sequence only, stripped of newlines
				f.write(fasta_record)
		


in_file_handle.close()
out_file_handle.close()

	