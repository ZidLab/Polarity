import sys
import os
import fnmatch
def features(destdir, input1path, input2dir):
    file1 = open(input1path,'r')
    names = input1path.split('/')
    names1 = names[len(names) - 1]
    filename = str(destdir) + str(names1)
    output = open(filename, 'w')
    file2s = []
    input2folder = os.listdir(input2dir)
    ext = '*.txt'
    for doc in input2folder:
        if fnmatch.fnmatch(doc, ext):
            docpath = str(input2dir) + str(doc)
            file2s.append(docpath)
    file2list = sorted(file2s, key=lambda x: int(x.split('D')[1].split('.')[0]))
    AlignCountTotal = 0.0
	
    print('begin')
    Chrcheck = 0
    for line in file1:
        temp = line.rstrip()
    	temp1 = temp.split('\t') 
    	Chr = int(temp1[1])
	if Chrcheck != 0:
		if Chrcheck != Chr:
			print('next chromosome' + str(Chr))
			Chrcheck = Chr
	else:
		print('first chromosome')
		Chrcheck = Chr
    	Posit = int(temp1[3])
    	Posit = Posit + 1
    	Orien = temp1[0]
    	count = temp1[2]
    	AlignCount = 0
    
    	#Seq = temp1[5]
    	if Chr == 1 :
            file2 = open(file2list[0],'r')
    	elif Chr == 2 :
            file2 = open(file2list[1],'r')
    	elif Chr == 3 :
            file2 = open(file2list[2],'r')
    	elif Chr == 4 :
            file2 = open(file2list[3],'r')
    	elif Chr == 5 :
            file2 = open(file2list[4],'r')
    	elif Chr == 6 :
            file2 = open(file2list[5],'r')
    	elif Chr == 7 :
            file2 = open(file2list[6],'r')
    	elif Chr == 8 :
            file2 = open(file2list[7],'r')
    	elif Chr == 9 :
            file2 = open(file2list[8],'r')
    	elif Chr == 10 :
            file2 = open(file2list[9],'r')
    	elif Chr == 11 :
            file2 = open(file2list[10],'r')
    	elif Chr == 12 :
            file2 = open(file2list[11],'r')
    	elif Chr == 13 :
            file2 = open(file2list[12],'r')
    	elif Chr == 14 :
            file2 = open(file2list[13],'r')
    	elif Chr == 15 :
            file2 = open(file2list[14],'r')
    	elif Chr == 16 :
            file2 = open(file2list[15],'r')
    	elif Chr == 17 :
            file2 = open(file2list[16],'r')    
    
    	for line2 in file2 :
            temp3 = line2.rstrip()
            temp2 = temp3.split('\t')
            Chr2 = int(temp2[1])
            if Chr == Chr2: #checking if the same chromosome
                if Orien == '-' :  # C orientation
                    name = temp2[0]
                    Orien2 = temp2[4]
                    if Orien2 == '-' :  #also C orientation
                        tStrt2 = int(temp2[3])
                    	tEnd2 = int(temp2[2])
                    	Strt2 = tStrt2+16
                    	End2 = tEnd2+13
                                    
                    	if End2 < Posit < Strt2 :
                            AlignCount = AlignCount + 1
                            output.write(name+'\t'+str(tStrt2)+'\t'+str(tEnd2)+'\t'+str(Posit)+'\t'+count+'\n') #sense count
                
                  
            	else :          # W orientation
                    name = temp2[0]
                    Orien2 = temp2[4]
                    if Orien2 == '+' :  #also W orientation
                        tStrt2 = int(temp2[2])
                   	tEnd2 = int(temp2[3])
                    	Strt2 = tStrt2-16
                    	End2 = tEnd2-13

                    	if Strt2 < Posit < End2 :
                            AlignCount = AlignCount + 1
                            output.write(name+'\t'+str(tStrt2)+'\t'+str(tEnd2)+'\t'+str(Posit)+ '\t'+count+'\n') #sense count
                	    AlignCountTotal = AlignCountTotal + float(count)
        file2.close()
    print AlignCountTotal
    print filename
    output.close()
    return filename
  
features(sys.argv[1], sys.argv[2], sys.argv[3])
