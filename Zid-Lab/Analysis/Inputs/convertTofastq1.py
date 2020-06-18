import sys
def convertTofastq1(path):
        file = open(path)
        names = path.split('.')
        filename = str(names[0]) + str('.fastq')
        print(filename)
        output = open(filename, 'w')
        n = 1
	inputs = []
	for line in file:
		if n%4 == 0:
                        inputs.append(line.rstrip())
                	#print(inputs)
			for i in range(4):
                        	output.write(inputs[i] + '\n')
			inputs = []
		else:
			if n%4 != 3:
				#print(inputs)
				inputs.append(line.rstrip())
			else:
				inputs.append('+')
		n += 1
	print(filename)
        output.close()
        return

print(sys.argv[1])
convertTofastq1(sys.argv[1])

