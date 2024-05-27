#### changes the names of reads to the file name to incorporate the depth and station of the read### Matt Hays summer 2019
#### Clara makes it so that it works for all the files in a folder
#### CD into the folder with the files

import glob
import os
path = r'./'

for metagenomicfilename in glob.glob(os.path.join(path, '*.fa')):
	file_name = metagenomicfilename
	metagenomicfilename = metagenomicfilename.replace('./','')
	metagenomicfilename = metagenomicfilename.replace('.fa','')
	outfilename = metagenomicfilename + '_name.fasta' 
	outfile = open( outfilename, 'w')

	#file_name_to_add = metagenomicfilename
	with open(file_name, 'r') as infile:
		for line in infile:
			if line.startswith('>'):
				splitline =  line.split(' ')
				readname = splitline[0]
				readname = readname.replace(":","_")
				gettingpair = splitline[1]
				gettingpair_split = gettingpair.split(":")
				pair = gettingpair_split[0]
				purereadname =  readname.split('>')
				outfile.write('>' + metagenomicfilename + '-' + purereadname[1] +'.'+pair+'\n')
			else:
				outfile.write(line)
	outfile.close()
			
			
		
















