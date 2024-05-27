infile=open('final.contigs.len1000.fa', 'r')
outfile=open('ETNP_1000m_final.contigs.len1000_name.txt', 'w')
for line in infile:
	if line.startswith('>'):
		line_lyst=line.split(' ')
		length_lyst=line_lyst[3].split('=')
		outfile.write(line_lyst[0]+'_L'+length_lyst[1])
	else:
		outfile.write(line)
outfile.close()
infile.close()