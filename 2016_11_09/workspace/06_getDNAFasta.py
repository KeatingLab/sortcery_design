#!/home/vxue/anaconda/bin/python

multiplexFile = "/home/vxue/data/experimental/SORTCERY/2016_11_09/multiplex.txt"
baseName = "/home/vxue/data/experimental/SORTCERY/2016_11_09/workspace/"
myMultiplex = open(multiplexFile,'r').readlines()


numIndices = int(len(myMultiplex)/13)
    

myNames = [(myMultiplex[i*13].strip()) for i in range(numIndices)]

import pandas as pd


bimPuma=['P','B']

for nameIdx in range(numIndices):
    allData = pd.DataFrame()

    for bg in [0,1]:
        foo = pd.read_csv(baseName+myNames[nameIdx]+"_dna/cell_counts_dnaframe_"+str(bg),delimiter='\s+',header=None,names=["CN_"+str(j) for j in range(12)])
        mySeqs = [seq.strip() for seq in open(baseName+myNames[nameIdx]+"_dna/unique_seqs_dnaframe_"+str(bg),'r').readlines()]
        foo['seq']=mySeqs
        foo['CN_tot']=foo.sum(axis=1)
        foo = foo[(foo.CN_tot>10) & ~foo.seq.str.contains("N")] #Sequences which contained N were thrown away
        foo['bg']=bimPuma[bg]

        if(bg==0):
            allData = foo
        else:
            allData = allData.append(foo)

    fastaOutDir = baseName+myNames[nameIdx]+"_dna/diversity.fasta"
    txtOutDir = baseName+myNames[nameIdx]+"_dna/diversity.txt"

    with open(fastaOutDir,'w') as writeStream, open(txtOutDir,'w') as txtWriteStream:
        for index,dnaSeq in allData.sort_values(by='CN_tot',ascending=False).reset_index(drop=True).iterrows():
            writeStream.write(">")
            writeStream.write(str(index))
            writeStream.write(";size=")
            writeStream.write(str(dnaSeq[["CN_"+str(j) for j in range(12)]].sum()))
            writeStream.write("\n")
            writeStream.write(dnaSeq.seq)
            writeStream.write("\n")

            txtWriteStream.write(str(dnaSeq[["CN_"+str(j) for j in range(12)]].sum()))
            txtWriteStream.write(" ")
            txtWriteStream.write(dnaSeq.seq)
            txtWriteStream.write("\n")
