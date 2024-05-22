normalize_dict = {}
outfile =open('Name_normalized_zeroed.txt', 'w')
gene_length = 1700 ######################################          USER CHANGES THE LENGTH FOR EACH GENE        ##################################
Sample_lyst = []
taxon_lyst = []
with open('Normalization_factor_file.csv', 'r') as infile:
	for line in infile:
		line_list = line.split(",")
		line_list[2] = line_list[2].replace('\n','')
		normalize_dict[line_list[0]] = line_list[2]
with open('Name_for_post_processing.txt', 'r') as infile:
	for line in infile:
		if line.startswith('name'):
			pass
		else:
			line_list = line.split(',')
			norm_factor = normalize_dict[line_list[0]]
			normalized_reads =float(line_list[2])*float(norm_factor)*100/gene_length
			line = line.replace('\n','')
			outfile.write(line + ',' + str(normalized_reads) + '\n')
			taxon = line_list[1]
			sample = line_list[0]
			if taxon in taxon_lyst:
				pass
			else:
				taxon_lyst.append(taxon)
			if sample in Sample_lyst:
				pass
			else:
				Sample_lyst.append(sample)

	for item in taxon_lyst:
		#print item
		y = []
		Difference=[]
		with open('Name_for_post_processing.txt', 'r') as infile:
			for line in infile:
				if line.startswith('name'):
					pass
				else:
					line_list = line.split(',')
					if line_list[1] == item:
						y.append(line_list[0])
					else: pass
			Difference = [t for t in Sample_lyst if t not in y]
			if len(Difference) > 0:
				for i in Difference:
					outfile.write(i+","+item+",0,0"+"\n")
			else:pass
outfile.close()
