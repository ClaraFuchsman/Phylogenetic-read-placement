
count = 0
outfile = open('Name.jplace', 'w')
with open('epa_result.jplace', 'r') as f:
	for line in f:
		if line.startswith('{\n'):
			outfile.write(line)
		elif line.startswith('  "tree'):
			outfile.write(line)
		elif line.startswith('  "placements"'):
			outfile.write(line)
		elif line.startswith('  [\n'):
			outfile.write(line + '    {"p": [')
			outfile.write('\n')
		elif line.startswith('      ['):
			list =  line.split(',')
			length =  list[4]
			length = length.replace(']' , '')
			length = float(length)
			if length <= 2.0:
				outfile.write(line)
				outfile.write(next(f))
				outfile.write(next(f))
				outfile.write(next(f))
				please = next(f)
				outfile.write(please)
				if please.startswith(','):
					outfile.write('    {"p": [')
			if length > 2.0:
				count = count +1
		elif line.startswith('  "metadata'):
			outfile.write(line)
		elif line.startswith('  "ver'):
			outfile.write(line)
		elif line.startswith('  "fields'):
			outfile.write(line)
		elif line.startswith('}'):
			outfile.write(line)
outfile.close()
print 'The number of reads above length 2 is...'
print count
		