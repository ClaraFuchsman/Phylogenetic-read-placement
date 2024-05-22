### This script takes the 'per query' file from the gappa result and counts the number of times the taxonomy occurs###
master_dict = {}
with open('per_query.tsv', 'r') as infile:
	for line in infile:
		line_list = line.split('	')
		sample_dragon = line_list[0]
		sample_list = sample_dragon.split('-')
		sample = sample_list[0]
		taxon = line_list[5].replace('\n', '')
		taxon_list = taxon.split(';')
		#want last member of lyst so that we don't double and triple count reads
		Taxonomy = taxon_list[-1]
		if sample + ',' + Taxonomy not in master_dict:
			master_dict[sample + ',' + Taxonomy] = 1
		else:
			 master_dict[sample + ',' + Taxonomy] += 1
import csv

with open('Name_for_post_processing.txt', 'w') as f:
    for key in master_dict.keys():
        f.write("%s,%s\n"%(key,master_dict[key]))









