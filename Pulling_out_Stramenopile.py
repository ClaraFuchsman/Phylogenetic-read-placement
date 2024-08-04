### This script takes the 'per query' file from the gappa result and pulls out the reads from a certain taxonomy###
outfile = open('Stramenopile_CB_reads.txt', 'w')
count = 0
master_lyst = []
with open('Protist_CB_details/per_query.tsv', 'r') as infile:
	for line in infile:
		line_list = line.split('	')
		read_name = line_list[0]
		taxon = line_list[5].replace('\n', '')
		taxon_list = taxon.split(';')
		#want last member of lyst so that we don't double and triple count reads
		Taxonomy = taxon_list[-1]
		if Taxonomy == str('Stramenopiles'):
			master_lyst.append(read_name)
		else:
			 pass
	#print master_lyst		 
fasta_filename = 'Protist_CB_reads.txt'			 
FASTA_file = open(fasta_filename).read()
FASTAs = FASTA_file.split('>')
for FASTA in FASTAs[1:]:
	FASTA_split = FASTA.split('\n')
	header = FASTA_split[0]
	if header in master_lyst:
		count = count + 1
		outfile.write(">" +FASTA)
	else: pass

print 'sequences in new file'
print count	
outfile.close()
infile.close()

		 
			 
			 
			 
			 
			 
			 


