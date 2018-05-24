#!/home/vxue/anaconda/bin/python

import pandas as pd
import os 


numBg = 2;
dirName = "/home/vxue/data/experimental/SORTCERY/2016_11_09/workspace/"
filePrefix = dirName+"seqframe/summary/seqframe_"
multiplexFile = "/home/vxue/data/experimental/SORTCERY/2016_11_09/multiplex.txt"
numBarcodes1 = 24;
numBarcodes2 = 2;

def main():
    myDataAssignments = getDataSplit(multiplexFile)
    for bg in range(numBg):
        for experiment in myDataAssignments:
            
            myExperiment = None
            
            for gate in myDataAssignments[experiment].keys():
                gateComponents = myDataAssignments[experiment][gate]
                
                sumGateCounts = None
                
                initialized=False
                for pair in gateComponents:
                    barcode1= pair[0]-1
                    barcode2 = pair[1]-1
                    barcode=(barcode2*24+barcode1)
                    
                    fileName = filePrefix+str(bg)+"_code_"+str(barcode)+".smry"
                    
                    print experiment,gate,"-",barcode1+1,barcode2+1,barcode
                    foo = pd.read_csv(fileName,
                    delimiter='\s+',
                    header=None,names=['bc'+str(barcode),'seq'])
                    
                    if(not initialized):
                        sumGateCounts = foo
                        initialized=True
                    else:
                        sumGateCounts = pd.merge(sumGateCounts,foo,on='seq',how='outer')
                
                sumGateCounts['CN_'+str(gate)]=sumGateCounts.sum(axis=1)#.info()
                
                if(myExperiment is None):
                    myExperiment = sumGateCounts
                else:
                    myExperiment = pd.merge(myExperiment,sumGateCounts,on='seq',how='outer')
                        
            myExperiment.fillna(0,inplace=True)
            #print myExperiment[['CN_'+str(each+1) for each in range(12)]].head()
            
            print experiment
            
            if not os.path.exists(dirName+experiment):
                os.makedirs(dirName+experiment)
            
            #Special case for Naieve
            numGates=12
            if (experiment=="Naive_Pool"):
                numGates=1
                
            
            with open(dirName+experiment+"/cell_counts_seqframe_"+str(bg),'w') as countStream, open(dirName+experiment+"/unique_seqs_seqframe_"+str(bg),'w') as seqStream:
                for index, each in myExperiment.iterrows():
                    seqStream.write(each.seq)
                    seqStream.write("\n")
                    for i in range(numGates):
                        countStream.write(str(int(each["CN_"+str(i+1)])) + " ")
                    countStream.write("\n")
            myExperiment.to_pickle(dirName+experiment+"/rawCounts_"+ str(bg) + ".pickle")

    
        
def getDataSplit(multiplexFile):
    
    mySamples = dict()
    with open(multiplexFile,'r') as readStream:
        index = 0
        
        myBarcodeLookup = None;
        
        for each in readStream:
            
            if index%13==0:                 #If line is a sample name, make a new dict
                myBarcodeLookup = dict() 
                mySamples[each.strip()]=myBarcodeLookup
            else:
                myBarcodeLookup[index%13]=[]
                    
            barcode1 = ""
            barcode2 = ""
            
    
    
            # For each pair
            index2=0
            for num in each.split():
                 
                if (index2==0): #skip gate index
                    index2+=1
                    continue
                    
                else:
    
                    if((index2+1)%2==0):
                        barcode1=num
                    else:
                        barcode2=num
                        myBarcodeLookup[index%13].append((int(barcode1),int(barcode2)))
                index2+=1        
     
            index+=1
    return mySamples

main()
