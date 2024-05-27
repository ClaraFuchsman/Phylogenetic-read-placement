infile=open('Name_taxon_2023.txt', 'r')
count_tab=0
count_line=0
for line in infile:
	count_line=count_line+1
	if '\t' in line:
		count_tab=count_tab+1
		splitline = line.split('\t')
		value = len(splitline)
		if len(splitline)==3:
			print line
		else: pass
	else:
		print "Here I am"
		print line
