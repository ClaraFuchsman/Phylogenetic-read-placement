import sys
from Bio import SeqIO
from Bio.Blast import NCBIXML
#### Clara makes it so that it works for all the files in a folder
#### CD into the folder with the files

import glob
import os
path = r'./'
outfile = open( 'Name_read_counts.txt', 'w')

for metagenomicfilename in glob.glob(os.path.join(path, '*.fasta')):
	file_handle = open(metagenomicfilename, 'rU')
	filename=metagenomicfilename.replace('./','')
	blastRecords = SeqIO.parse(file_handle, 'fasta')
	counter = 0

	for rec in blastRecords:
		counter += 1
	
	print filename + ' '+ str(counter)+'\n'
	outfile.write(filename + ' '+ str(counter)+'\n')
outfile.close()