########IMPORTANT######## DO NOT USE A .fasta FILE AS YOUR REFERENCE FILE###################################
taxon_search = []

with open('list_of_taxons.txt', 'r') as f:
	for line in f:
		line = line.replace('\n','')
		taxon_search.append(line)
		for i in taxon_search:
			i_ = i.replace(" ","_")
			filename = "%s.txt" %i_

for i in taxon_search:
	print i
	i_ = i.replace(" ","_")
	filename = "%s.txt" %i_
	outfile = open(filename, 'a+')
	with open('Monirussaman_taxon.txt','r') as infile:
		for line in infile:
			if i in line:
				line_list = line.split('\t')
				outfile.write(line_list[0] + '\n')
				print line_list[0]
			else:pass
	outfile.close()
	############################
	# enter the name of the original .fasta file used to make the tree:
	fasta_filename = 'PolB_slim_more.txt'

	# enter the name of the text file containing the sequences to be kept - one per line
	print filename
	seqs_to_be_kept = filename
	# must have unix linebreaks

	# enter the name of the new .fasta file
	outfile_name = '%s_sequence.fa' %i_
	###############################

	FASTA_file = open(fasta_filename).read()
	FASTAs = FASTA_file.split('>')

	keep_list = []
	textfile = open(seqs_to_be_kept)
	for line in textfile:
		if len(line) < 2: pass
		else:
			name = line.strip('\n')
			name = name.strip('>')
			keep_list.append(name)
	#print keep_list

	outfile2 = open(outfile_name, 'a+')

	count = 0
	for FASTA in FASTAs[1:]:
		FASTA_split = FASTA.split('\n')
		header = FASTA_split[0]
		header = header.split(' ')
		header = header[0]
		header = header.split('\t')
		header = header[0]
		header = header.split('.')
		header = header[0]

		#print header
		if header in keep_list:
			outfile2.write('>')
			outfile2.write(FASTA)
			count = count + 1
	#print count,
	#print 'output sequences'
	#print len(keep_list),
	#print 'sequences in to keep file'

	textfile.close()				
	outfile2.close()		

import glob
import os
path = r'./'
#from __future__ import print_function
from Bio import pairwise2 as pw2
from Bio import SeqIO
master_percent = open('similarity_averages.txt', 'w')
for ref_file in glob.glob(os.path.join(path, '*.fa')):
	dictionary={}
	count = 0
	percent_similarity_sum = 0
	out_name = str(ref_file)
	out_name = out_name.replace('sequence.fa','sim.txt')
	out_name = out_name.replace('./','')
	taxon_sim_file = open(out_name, 'w')
	for record in SeqIO.parse(ref_file,"fasta"):
		dictionary[record.id] = record.seq
			
	for record in SeqIO.parse(ref_file,"fasta"):
		#count = 0
		#percent_similarity_sum = 0
		#print(out_name)
		for i in dictionary:
			seq_length = min(len(record.seq), len(dictionary[i]))
			global_align = pw2.align.globalxx(record.seq, dictionary[i])
			matches = global_align[0][2]
			percent_match = (matches / seq_length) * 100
			percent_similarity_sum += percent_match
			taxon_sim_file.write(record.id + ',' + i + ',' + str(percent_match) + '\n')
			count += 1
	taxon_sim_file.close()
	average_similarity = percent_similarity_sum / count
	master_percent.write(out_name + ',' + str(average_similarity) + '\n')
			
	print count





















