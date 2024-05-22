#pretty much all this script does is add a > so it can run through epa-ng after going through papara. for reads file
outfile = open('Name_READS_spacedash.txt', 'w')
with open('Name_READS_joinedReadsForPlacement.txt', 'r') as infile:
	next(infile)
	for line in infile:
		if ' -' in line:
			line_list = line.split(' -')
			line_list[0] = line_list[0].replace(' ', '')
			outfile.write('>' + line_list[0] + '\n')
			outfile.write('-' + line_list[1])
		elif ' A' in line:
			line_list = line.split(' A')
			line_list[0] = line_list[0].replace(' ', '')
			outfile.write('>' + line_list[0] + '\n')
			outfile.write('A' + line_list[1])
		elif ' B' in line:
			line_list = line.split(' B')
			line_list[0] = line_list[0].replace(' ', '')
			outfile.write('>' + line_list[0] + '\n')
			outfile.write('B' + line_list[1])
		elif ' C' in line:
			line_list = line.split(' C')
			line_list[0] = line_list[0].replace(' ', '')
			outfile.write('>' + line_list[0] + '\n')
			outfile.write('C' + line_list[1])
		elif ' D' in line:
			line_list = line.split(' D')
			line_list[0] = line_list[0].replace(' ', '')
			outfile.write('>' + line_list[0] + '\n')
			outfile.write('D' + line_list[1])
		elif ' E' in line:
			line_list = line.split(' E')
			line_list[0] = line_list[0].replace(' ', '')
			outfile.write('>' + line_list[0] + '\n')
			outfile.write('E' + line_list[1])
		elif ' F' in line:
			line_list = line.split(' F')
			line_list[0] = line_list[0].replace(' ', '')
			outfile.write('>' + line_list[0] + '\n')
			outfile.write('F' + line_list[1])
		elif ' G' in line:
			line_list = line.split(' G')
			line_list[0] = line_list[0].replace(' ', '')
			outfile.write('>' + line_list[0] + '\n')
			outfile.write('G' + line_list[1])
		elif ' H' in line:
			line_list = line.split(' H')
			line_list[0] = line_list[0].replace(' ', '')
			outfile.write('>' + line_list[0] + '\n')
			outfile.write('H' + line_list[1])
		elif ' I' in line:
			line_list = line.split(' I')
			line_list[0] = line_list[0].replace(' ', '')
			outfile.write('>' + line_list[0] + '\n')
			outfile.write('I' + line_list[1])
		elif ' K' in line:
			line_list = line.split(' K')
			line_list[0] = line_list[0].replace(' ', '')
			outfile.write('>' + line_list[0] + '\n')
			outfile.write('K' + line_list[1])
		elif ' L' in line:
			line_list = line.split(' L')
			line_list[0] = line_list[0].replace(' ', '')
			outfile.write('>' + line_list[0] + '\n')
			outfile.write('L' + line_list[1])
		elif ' M' in line:
			line_list = line.split(' M')
			line_list[0] = line_list[0].replace(' ', '')
			outfile.write('>' + line_list[0] + '\n')
			outfile.write('M' + line_list[1])
		elif ' N' in line:
			line_list = line.split(' N')
			line_list[0] = line_list[0].replace(' ', '')
			outfile.write('>' + line_list[0] + '\n')
			outfile.write('N' + line_list[1])
		elif ' P' in line:
			line_list = line.split(' P')
			line_list[0] = line_list[0].replace(' ', '')
			outfile.write('>' + line_list[0] + '\n')
			outfile.write('P' + line_list[1])
		elif ' Q' in line:
			line_list = line.split(' Q')
			line_list[0] = line_list[0].replace(' ', '')
			outfile.write('>' + line_list[0] + '\n')
			outfile.write('Q' + line_list[1])
		elif ' R' in line:
			line_list = line.split(' R')
			line_list[0] = line_list[0].replace(' ', '')
			outfile.write('>' + line_list[0] + '\n')
			outfile.write('R' + line_list[1])
		elif ' S' in line:
			line_list = line.split(' S')
			line_list[0] = line_list[0].replace(' ', '')
			outfile.write('>' + line_list[0] + '\n')
			outfile.write('S' + line_list[1])
		elif ' T' in line:
			line_list = line.split(' T')
			line_list[0] = line_list[0].replace(' ', '')
			outfile.write('>' + line_list[0] + '\n')
			outfile.write('T' + line_list[1])
		elif ' V' in line:
			line_list = line.split(' V')
			line_list[0] = line_list[0].replace(' ', '')
			outfile.write('>' + line_list[0] + '\n')
			outfile.write('V' + line_list[1])
		elif ' W' in line:
			line_list = line.split(' W')
			line_list[0] = line_list[0].replace(' ', '')
			outfile.write('>' + line_list[0] + '\n')
			outfile.write('W' + line_list[1])
		elif ' X' in line:
			line_list = line.split(' X')
			line_list[0] = line_list[0].replace(' ', '')
			outfile.write('>' + line_list[0] + '\n')
			outfile.write('X' + line_list[1])
		elif ' Y' in line:
			line_list = line.split(' Y')
			line_list[0] = line_list[0].replace(' ', '')
			outfile.write('>' + line_list[0] + '\n')
			outfile.write('Y' + line_list[1])
		elif ' Z' in line:
			line_list = line.split(' Z')
			line_list[0] = line_list[0].replace(' ', '')
			outfile.write('>' + line_list[0] + '\n')
			outfile.write('Z' + line_list[1])
outfile.close()








































