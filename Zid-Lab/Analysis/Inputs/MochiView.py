import os
def MochiView(destDir, inputpath)
    names = inputpath.split('.')
    filename = str(destDir) + str(names[len(names) - 2]) + str('.txt')
    file1 = open(inputpath,'r')
    output = open(filename, 'w')
    output.write('SEQ_NAME\tSTART\tEND\tSTRAND\tCOUNT\n')
    #temp7 = file1.readline()
    norm = 2.65 # Normalized to total alignable reads/3 million

    for line in file1:
        temp = line.rstrip()
        temp1 = temp.split('\t') 
        name = "Chr"+temp1[1]
       
        PositTemp = int(temp1[3])
        Orien = temp1[0]
    	if Orien == '+' :
        	count = (float(temp1[2]))/norm
        	Start = str(PositTemp)
        
    	else :
        	tempc = (float(temp1[2]))/norm
        	count = tempc*(-1)
        	Start = str(PositTemp)
        
    output.write(name+'\t'+Start+'\t'+Start+'\t'+Orien+'\t'+str(count)+ '\n')  #antisense count
    output.close()
    file1.close()
    return filename
