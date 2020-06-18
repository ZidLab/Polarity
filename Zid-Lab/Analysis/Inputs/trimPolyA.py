#!/bin/python
import sys
import os
import fnmatch
def trimPolyA(filename,inputpath):
	file1 = open(inputpath)
	print(filename)
	output = open(filename, 'w+')
	count = 0
	count1 = 0
	count2 = 0
	num = 0
	read = 0
	for line in file1 :
		num+=1
		if num ==2:  #The file has sequences on the 2nd line, every 4 lines
			temp2 = line.rstrip() #?
			temp26 = temp2[0:25] #Take the first 26 bps
			tempminus1 = temp26.rstrip('A') #Strip all the A's
			tempminus2 = tempminus1.rstrip('N') #Strip the next base (in case an A was misread)
			tempminusA = tempminus2.rstrip('A') #Strip more A's
			Length = len(tempminusA)
			tempTest = tempminusA
			tempTest1 = tempTest[0:-1]
			tempTest2 = tempTest1.rstrip('A')
			LenChange = Length - len(tempTest2) #I'm double checking for extra runs of A's, not sure if its needed.
			if LenChange > 6 :
				count2 = count2 + 1
				read = 7
			else:
				if Length > 17:
					if Length > 23:
						Sequence = tempminusA[0:23] #Only use the first 23bps for alignment 0:22
						Length = len(Sequence)      #More likely to have misreads later in the sequence
					else:
						Sequence = tempminusA
					acount = Sequence.count('A')
					temp4 = Length - 5
					if acount < temp4 :
						count = count + 1
						read = 0
					else:
						read = 7
				else:
					read = 7

		elif num == 4: #This is to get the Quality of the sequence read (num%4 == 2)
			count1+=1
			num = 0 #n don't need this
			QualityT = line.rstrip()
			Quality = QualityT[0:Length]
			if read == 7:
				read = 0
			else:
				output.write('@r'+str(count)+'/1.'+str(Length)+'\n'+Sequence+'\n'+'+\n'+Quality+'\n')


	output.close()
	print("trimmed")
	#Takes txt Solexa Sequence results and outputs them in fastq for
	#Took the first 26bp, removed the 3'AAAAs
	return
trimPolyA(sys.argv[1],sys.argv[2])
