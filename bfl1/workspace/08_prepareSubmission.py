import os 

numBg = 1;
dirName = "/home/vxue/data/SORTCERY_PUBLICATION/bfl1/workspace/"
filePrefix = dirName+"seqframe/summary/seqframe_"
multiplexFile = "/home/vxue/data/SORTCERY_PUBLICATION/bfl1/multiplex.txt"
numBarcodes1 = 24;
numBarcodes2 = 2;

def main():
    myDataAssignments = getDataSplit(multiplexFile)
    print("experiment", "unique_barcode_identifier", "sortcery_gate", "barcode1", "barcode2")
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
                    print( experiment, str(barcode) , str(gate), barcode1, barcode2)
		    
            
            #Special case for Naieve
            numGates=12
            if (experiment=="Naive_Pool"):
                numGates=1
                
    
        
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
