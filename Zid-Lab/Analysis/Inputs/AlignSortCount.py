import sys
def AlignSortCount(inputpath, destdir, mappath):
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
	print('begin')
	file1 = open(mappath,'r')
	names = inputpath.split('/')
	filename = str(destdir) + names[len(names) - 1]
	output = open(filename, 'w')
	count2 = 0
    	count3 =0
    	temp = file1.readline()
    	line1 = temp.rstrip() 
	temp1 = line1.split('\t')
	ChrTemp2 = temp1[2]	
	if ChrTemp2 == "MT":
        	Chr = 17
	else:
        	Chr = roman_to_arabic(ChrTemp2)

    	Name = temp1[0]
    	NameLen = Name[-2:]
    	Orien = temp1[1]
    	if Orien == '-':
        	Posit = (int(temp1[3]))+(int(NameLen))  #didn't subtract 1 adjusting because it is one off
    	else :
        	Posit = (int(temp1[3]))+1 #adjusting because it seems to be one off on the +Strand


    	dct_lst = [{"Chr":Chr,"Posit":Posit,"Orien": Orien}]#"Seq": 'a'}]
    	for line in file1 :
        	count3 +=1
        	line1 = line.rstrip() 
        	temp1 = line1.split('\t')
        	ChrTemp2 = temp1[2]
        	if ChrTemp2 == "MT":
            		Chr = 17
        	else:
            		Chr = roman_to_arabic(ChrTemp2)
        	Name = temp1[0]
        	NameLen = Name[-2:]
        	Orien = temp1[1]
        	if Orien == '-':
            		Posit = (int(temp1[3]))+(int(NameLen))
        	elif Orien =='+' :
            		Posit = (int(temp1[3]))+1
        	temp_lst = [{"Chr":Chr,"Posit":Posit,"Orien": Orien}]#"Seq": 'a'}]#temp1[4]
        	dct_lst.extend(temp_lst)

    	from operator import itemgetter
    	dct_lst.sort(key=itemgetter("Posit"))

    	from operator import itemgetter
    	dct_lst.sort(key=itemgetter("Chr"))

    	count = 1
    	Length = len(dct_lst)
    	for i in range(len(dct_lst)) :
        	Chr = (dct_lst[i]["Chr"])
        	Posit = (dct_lst[i]["Posit"])
        	Orien = (dct_lst[i]["Orien"])
        	#Seq = (dct_lst[i]["Seq"])
        	if i == (Length-1):
            		output.write(Orien+'\t'+str(Chr)+'\t'+str(count)+'\t'+str(Posit)+'\t'+'\n')#Seq+'\n')
        	else:
            		OrienT2 = (dct_lst[i+1]["Orien"])
            		PositT2 = (dct_lst[i+1]["Posit"])
            		if Orien == OrienT2 :
        
                		if Posit == PositT2 :
                    			count = count + 1
                		else :
                    			count2 = count + count2
                    			output.write(Orien+'\t'+str(Chr)+'\t'+str(count)+'\t'+str(Posit)+'\t'+'\n')#Seq+'\n') #sense count
                    			count = 1
            		else :
                		count2 = count + count2
                		output.write(Orien+'\t'+str(Chr)+'\t'+str(count)+'\t'+str(Posit)+'\t'+'\n')#Seq+'\n') #sense count
                		count = 1
         
    	output.close()  
	file1.close()
    	print count2	
	print filename
	return
	
	
AlignSortCount(sys.argv[1], sys.argv[2], sys.argv[3])
	
	
	
	
	
	
	
	
	
	

	

	
		

	
	
		
	
	
	
	
	
	
	
	

