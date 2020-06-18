import sys
import pandas as pd
import numpy as np
import os
def calculatePolarity(i, l, densityati, density):
    normDist = (2*i - (l + 1))/(l - 1)   
    polarity = densityati*normDist/density
    
    return polarity

def PolarityPerGene(inputpath, destdir):
    df = pd.read_csv(inputpath)
    df = df.iloc[:,-8:]
    df = df[df.RelativePosition >= 15]
    names = inputpath.split('.')
    names1 = names[0].split('/')
    #this is the final output file
    nameofcsv = destdir + names1[len(names1) - 1] + str("_polarity&gene.csv")
    print(nameofcsv)
    
    len_dict = dict(zip(df['Name'],df['Length']))
    den_dict = df.groupby('Name')['Counts'].agg('sum').to_dict()
    
    Polaritylist = []
    for i in range (len(df)):
        polarityofi = calculatePolarity(df.iloc[i,6], df.iloc[i,5], df.iloc[i,7], den_dict.get(df.iloc[i,0]))
        Polaritylist.append(polarityofi)
    df['Polarity'] = Polaritylist
    
    #calculate polarity per gene
    pol_dict = df.groupby('Name')['Polarity'].agg('sum').to_dict()
    
    Lengthlist = []
    Totalcountlist = []
    names = []
    genepol = []
    for name in pol_dict.keys():
        if name in len_dict.keys() and name in den_dict.keys():
            Lengthlist.append(len_dict.get(name))
            Totalcountlist.append(den_dict.get(name))
            names.append(name)
            genepol.append(pol_dict.get(name))
        
    d = {'Name':names,'Length':Lengthlist,'Counts':Totalcountlist,'PolarityPerGene': genepol}
    dfcsv = pd.DataFrame(data=d)
    dfcsv = dfcsv.dropna(how='any',axis=0)
    dfcsv.to_csv(nameofcsv)

inputlist = []
inputfolder = sys.argv[1]
ext = "fastq"
for root, dirs, files in os.walk(inputfolder):
    for file in files:
        if file.split('.')[-1] == ext:
            inputlist.append(inputfolder+file)
print('\n Number of input files: ' + str(len(inputlist))) 
for i in range(0, len(inputlist)):
    PolarityPerGene(inputlist[i], sys.argv[2])
