import sys
def countRibo(destdir, inputpath):
	numerals = [{'letter': 'M', 'value': 1000}, {'letter': 'D', 'value': 500}, {'letter': 'C', 'value': 100}, {'letter': 'L', 'value': 50}, {'letter': 'X', 'value': 10}, {'letter': 'V', 'value': 5}, {'letter': 'I', 'value': 1}]
	def roman_to_arabic(number):
		index_by_letter = {}
		for index in xrange(len(numerals)):
			index_by_letter[numerals[index]['letter']] = index
		result = 0
		previous_value = None
		for letter in reversed(number):
			index = index_by_letter[letter]
			value = numerals[index]['value']
			if (previous_value is None) or (previous_value <= value):
				result += value
			else:
				result -= value
			previous_value = value
		return result

	file1 = open(inputpath,'r')
	names = inputpath.split('/')
	names1 = names[len(names) - 1].split('.s')
	#temp = inputpath.split('.')[0] + str('-temp.csv')
	filename = str(destdir) + names1[len(names1) - 2] + str('-aligned.txt')
	#output = open(temp, 'w')
	output = open(filename, 'w')
	print(filename)
	for line in file1:
		if not line.startswith('@'):
			strs = line.split('\t')
			inputs = []
			if str(strs[1]) == '16':
				inputs.append('-')
			else:
				inputs.append('+')
			chromnum = strs[2][3:]
        		if chromnum == "MT":
                		Chr = 17
        		else:
                		Chr = roman_to_arabic(chromnum)
			inputs.append(str(Chr))
			inputs.append('1')
			inputs.append(strs[3])
			#print(inputs)
			for i in range(len(inputs)):
				output.write(inputs[i])
				output.write('\t')
			output.write('\n')	
	return			
countRibo(sys.argv[1], sys.argv[2])
