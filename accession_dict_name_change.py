###used to take the Prokka output table and match all of the gene locus tags with the accession number###
import re

id_dict = {}

Identifier ='ETNP_1000m'
fp = open('ETNP_1000m.tbl', 'r')

for line in fp:
    if line.startswith('>'):
        #print line
        a_list = line.split('>')
        accession = a_list[1]
        accession = accession.replace('\n', '')
        accession = accession.replace('Feature ', '')
        #print (accession)
    elif 'locus_tag' in line:
        l_list = line.split('locus_tag')
        #print l_list
        locus_tag = l_list[1]
        locus_tag = locus_tag.replace('\n', '')
        locus_tag = locus_tag.replace('\t', '')
        #print locus_tag
        id_dict[locus_tag] = accession
    else:
        pass
        
fp.close()
#print (id_dict)
#######Changes the name of the .faa file of the prokka output#######
outfile = open('ETNP_1000m_prokka_name.faa', 'w')
with open('ETNP_1000m.faa', 'r+') as f:
	for line in f:
		if line.startswith('>'):
			list = line.split(' ')
			line = line.replace(' ','_')
			line = line.replace('>','')
			a = list[0]
			a = a.replace('>','')
			b = id_dict[a]
			line = '>' + Identifier + '_' + b + '-' + line
			outfile.write(line)
			#print line
		else:
			outfile.write(line)
outfile.close()



















































