#### changes the names of reads to the file name to incorporate the depth and station of the read### Matt Hays summer 2019
#### Clara makes it so that it works for all the files in a folder
#### CD into the folder with the files

import glob
import os
path = r'./'

for metagenomicfilename in glob.glob(os.path.join(path, '*.fasta')):
	file_name = metagenomicfilename
	metagenomicfilename = metagenomicfilename.replace('./','')
	metagenomicfilename = metagenomicfilename.replace('.fasta','')
	outfilename = metagenomicfilename + '_name.fasta' 
	outfile = open( outfilename, 'w')

	file_name_to_add = metagenomicfilename.split('.')
	with open(file_name, 'r') as infile:
		for line in infile:
			if line.startswith('>'):
				splitline =  line.split(' ')
				readaccession = splitline[0]
				purereadaccession =  readaccession.split('>')
				outfile.write('>' + file_name_to_add[0] + '-' + purereadaccession[1]+'.'+file_name_to_add[2]+'\n')
			else:
				outfile.write(line)
	outfile.close()
			
			
		
















