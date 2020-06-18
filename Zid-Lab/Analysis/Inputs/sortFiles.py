import sys
def sortFile(destdir, inputpath):
    file1 = open(inputpath,'r')
    names = inputpath.split('/')
    names1 = names[len(names) - 1]
    filename = str(destdir) + str(names1) + str('.csv')
    output = open(filename, 'w')
    temp = file1.readline()
    line1 = temp.rstrip()
    temp1 = line1.split('\t')
    Name = temp1[0]
    ChrT = Name[2]
    Chr = '1'
    dct_lst = [{"name":temp1[0],"Start":temp1[1],"End":"7236","Posit":temp1[3],"Count":temp1[4]}]

    print('begin')
    #output.write('Name'+'\t'+'Start'+'\t'+'End'+'\t'+'Posit'+ '\t'+'Count'+'\n')

    for line in file1 :
        temp = line.rstrip()
        temp1 = temp.split('\t')

        temp_lst = [{"name":temp1[0],"Start":temp1[1],"End":temp1[2],"Posit":temp1[3],"Count":temp1[4]}]
        dct_lst.extend(temp_lst)
    file1.close()

    from operator import itemgetter
    dct_lst.sort(key=itemgetter("name"))


    for i in range(len(dct_lst)) :
        output.write(dct_lst[i]["name"]+','+dct_lst[i]["Start"]+','+dct_lst[i]["End"]+','+dct_lst[i]["Posit"]+','+dct_lst[i]["Count"]+'\n') #sense count

    output.close()
    print(filename)

sortFile(sys.argv[1], sys.argv[2])
