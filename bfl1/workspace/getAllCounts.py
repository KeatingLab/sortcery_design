#!/home/vxue/anaconda/bin/python

import pandas as pd


numBg = 2;
dirName = "/home/vxue/luther/SORTCERY/SORTCERY2_all/"
filePrefix = dirName+"seqframe/summary/seqframe_"
numBarcodes = 24;


for bg in range(numBg):
    myDataFrame = None;

    for barcode in range(numBarcodes):
        fileName = filePrefix+str(bg)+"_code_"+str(barcode)+".smry"
            
        foo = pd.read_csv(fileName,
                    delimiter='\s+',
                    header=None,names=['cn'+str(barcode),'seq'])
        
        
        
        if(barcode==0):
            myDataFrame = foo
        else:
            myDataFrame = pd.merge(myDataFrame,foo,on='seq',how='outer')
            
    myDataFrame.fillna(0,inplace=True)
    
    outputDF = myDataFrame[['seq']].copy();
    for each in range(12):
        outputDF['CN'+str(each)] = myDataFrame['cn'+str(each)]+myDataFrame['cn'+str(each+12)]
    

    
    with open(dirName+"cell_counts_seqframe_"+str(bg),'w') as countStream, open(dirName+"unique_seqs_seqframe_"+str(bg),'w') as seqStream:
        for index, each in outputDF.iterrows():
            seqStream.write(each.seq)
            seqStream.write("\n")
            for i in range(12):
                countStream.write(str(int(each["CN"+str(i)])) + " ")
            countStream.write("\n")