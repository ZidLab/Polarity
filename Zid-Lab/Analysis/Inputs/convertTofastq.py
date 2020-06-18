import sys
def convertTofastq(path):
	file = open(path)
	names = path.split('.')
	filename = str(names[0]) + str('.fastq')
	print(filename)
	output = open(filename, 'w')
	for line in file:
		strs = line.split('/1:')
		inputs = []
		inputs.append(strs[0] + '/1')
		strs1 = strs[1].split(':')
		for j in range(len(strs1)):
			inputs.append(strs1[j])
		#print(inputs)
		for i in range(4):
			if i == 2:
				output.write('+' + '\n')
			elif i == 3:
				output.write(inputs[i - 1])
			else:
				output.write(inputs[i] + '\n')
	print(filename)
	output.close()
	return 

print(sys.argv[1])
convertTofastq(sys.argv[1])	
