#!/home/vxue/anaconda/bin/python

from Bio import SeqIO,Seq
import multiprocessing
import sys 

'''
This code is a little difficult 
important parameters to update are 

myUniqueName

'''
myUniqueName = "2016_11_09"



globalNumMutations = 5   # How many mutations are allowed in the overlapping region
globalNumMaxOverlap = -1 # 0 for align to ref 
                         #-1 for test all sliding positions
globalMinNumOverlap = 70  # Minimum number of overlap before mapping to ref
globalReference=[]

#global globalReference
globalReference =["",""]
globalReference[0] = "GGAGGCGGTAGCGGAGGCGGAGGGTCGGCTAGCGGTACCGGATCCGGTGGCCAATGGGCGCGTGAAATTGGCGCCCAACTGCGTCGCATGGCGGATGATCTGAATGCCCAATATGAACGTCGTCGCCAGGAGGAACAACAAAGAGGCGGCCGCGAT"
globalReference[1] = "GGAGGCGGTAGCGGAGGCGGAGGGTCGGCTAGCGGTACCGGATCCGGTGGCCGTCCGGAAATTTGGATTGCGCAGGAACTGCGTCGTATTGGCGATGAATTTAATGCGTATTATGCGCGTCGCGTGTTTCTGAATAACTATCAGGGCGGCCGCGAT"


globalVarPos=[[],[]]
globalVarPos[0]=[54,57,60,63,66,69,75,81,87,90,96,99,102,105,111]
globalVarPos[1]=[54,57,60,63,66,69,75,81,87,90,96,99,102,105,111]

globalRangeToExtract=[48,114]

globalReferenceScaffold=["",""]

globalNCores = 16

writeQueue = None

def main(barcodeToProcess):
    
    dirName = "/home/vxue/data/experimental/SORTCERY/"+myUniqueName+"/workspace/"
    initializeScaffold()

    processFile(dirName+"barcode_",barcodeToProcess,dirName)
    
def initializeScaffold():
    global globalReferenceScaffold
    globalReferenceScaffold[0]=globalReference[0]
    globalReferenceScaffold[1]=globalReference[1]

    # Do you want the reference alignment to be 
    # ...AAA...DDD......GGG (Use snippet below)
    # AAABBBCCCDDDEEEFFFGGG (Use snippet above)
    # 

    #for bg in range(len(globalReference)):
    #    myRef = list(globalReference[bg])
    #    for aaPos in globalVarPos[bg]:
    #        myRef[aaPos]='.'
    #        myRef[aaPos+1]='.'
    #        myRef[aaPos+2]='.'
    #        
    #    globalReferenceScaffold[bg] = "".join(myRef)
    
def processFile(prefix,barcodeToProcess,dirName):

    records = SeqIO.parse(open(prefix+str(barcodeToProcess),"rU"), "fastq-sanger")


    manager = multiprocessing.Manager()
    global writeQueue
    writeQueue = manager.Queue()
    pool = multiprocessing.Pool(processes=globalNCores)

    #put listener to work first
    writeListener = pool.apply_async(writeToFile, (writeQueue,dirName,barcodeToProcess))
    
    
    myIt = pool.imap(func=getDNASeq,iterable=pairIterator(records),chunksize=1000)
    
    for each in myIt:
        each
    

    writeQueue.put('kill')
    pool.close()
    print("COMPLETE!")

def pairIterator(records):
    for pair1,pair2 in zip(records,records):
        yield (pair1,pair2)

def writeToFile(q, dirName,barcode):
    '''listens for messages on the q, writes to file. '''
    fileArray1 = [open(dirName+"/seqframe/seqframe_"+str(each)+"_code_"+barcode+"",'w') for each in range(len(globalReference))]
    fileArray2 = [open(dirName+"/dnaframe/dnaframe_"+str(each)+"_code_"+barcode+"",'w') for each in range(len(globalReference))]
    fileArray3 = [open(dirName+"/qualframe/qualframe_"+str(each)+"_code_"+barcode+"",'w') for each in range(len(globalReference))]

    while 1:
        res = q.get()
        if res == 'kill':
            break
        #refIndex,extractedAA,fullDNA,str(extractedQuality)
        f = fileArray1[res[0]]    
        f.write(str(res[1]) + '\n')
        f.flush()
        
        f = fileArray2[res[0]]    
        f.write(str(res[2]) + '\n')
        f.flush()
        
        f = fileArray3[res[0]]    
        f.write(str(res[3]) + '\n')
        f.flush()
        
    
    [each.close() for each in fileArray]
    [each.close() for each in fileArray2]
    [each.close() for each in fileArray3]
    
def passQuality(forwardRecord,reverseRecord):
    forward_Quality = forwardRecord.letter_annotations['phred_quality']
    reverse_Quality = reverseRecord.letter_annotations['phred_quality'][::-1]
    for i in range(40):
        if(forward_Quality[i]<20 or reverse_Quality[i] < 20):
            return False
    return True

def getAlignmentsForBestMatch(string1,string2):
    #Sliding window matcher
    
    extendedString1 = "."*len(string2)+string1+"."*len(string2);
    extendedString2 = string2
    counter=0
    shiftArray=[]
    #For all positions in the reference
    for i in range(len(string1)+len(string2)):
        counter = 0;
        numCompared=0;

        for j in range(len(extendedString2)):
            if(extendedString2[j]=='.' or extendedString1[j]=='.'):
                counter+=0
            elif(extendedString2[j]==extendedString1[j]):
                counter+=1;
                numCompared+=1
            else:
                numCompared+=1
                
        #print extendedString1
        #print extendedString2 , counter
        #print '----------'
        shiftArray.append((i-len(string2),counter,numCompared))    
        extendedString2 = "."+extendedString2;

    return sorted(shiftArray,key = lambda x:x[1],reverse=True)
    
def getAlignmentsForTailMatch(string1,string2,maxOverlap=-1,allowableMutations=0,verbose=False):
    #Sliding window matcher
    
    extendedString1 = string1+"."*len(string2); #Extend the first string to accomadate second   #AAAAA.........
    extendedString2 = "."*len(string1)+string2  #Move the seond string to the tail of the first #.....BBBBBB
    counter=0
    shiftArray=[]
    
    thisRange=0
    if(maxOverlap==-1):
        thisRange = len(string1)
    else:
        thisRange=maxOverlap
    
    
    #For all positions in the reference
    for i in range(thisRange+1):
        counter = 0;
        numCompared=0;

        aString = ""
        bString = ""
        
        for j in range(len(extendedString2)):
            if(extendedString2[j]=='.' or extendedString1[j]=='.'):
                counter+=0
            elif(extendedString2[j]==extendedString1[j]):
                counter+=1;
                numCompared+=1
                aString=aString+extendedString1[j]
                bString=bString+extendedString2[j]
            else:
                numCompared+=1
                aString=aString+extendedString1[j]
                bString=bString+extendedString2[j]
        
        
        if (verbose):
            #print aString,bString
            print(extendedString1)
            print(extendedString2 , scoreAlignment(aString,bString) , i)
            #print '----------'
        
        shiftArray.append((len(string1)-i,counter,numCompared,scoreAlignment(aString,
                                                                             bString,
                                                                             allowableMutations=allowableMutations)))    
        extendedString2 = extendedString2[1:];

    return sorted(shiftArray,key = lambda x:x[3],reverse=True)

def scoreAlignment(aString,bString,allowableMutations=0):
    #Takes in two aligned strings and scores them
    score=allowableMutations
    slack=allowableMutations
    
    for i in range(len(aString)-1,-1,-1):
        #print aString[i],aString[i]
        if(aString[i]==bString[i]):
            score+=1
        else:
            slack-=1
            if(slack==-1):
                return -1
    return score



def getAASeq(myInput):
    forwardRecord,reverseRecord=myInput
    if(passQuality(forwardRecord,reverseRecord)):
        refIndex,fullDNA,extractedAA = getDNAAndAA(forwardRecord,reverseRecord)
        
        expectedLength = (globalRangeToExtract[1]-globalRangeToExtract[0])/3
        if(expectedLength == len(extractedAA)):
            writeQueue.put((refIndex,extractedAA))
            return 1
    return 0

def getDNASeq(myInput):
    forwardRecord,reverseRecord=myInput
    if(passQuality(forwardRecord,reverseRecord)):
        refIndex,fullDNA,extractedAA,extractedQuality = getDNAAndAA(forwardRecord,reverseRecord)
        
        expectedLength = globalRangeToExtract[1]-globalRangeToExtract[0]
        if(expectedLength == len(fullDNA)):
            writeQueue.put((refIndex,extractedAA,fullDNA,str(extractedQuality)))

            return 1
    return 0

def getBestBG(myRead):
    bestBG = -1
    bestNumDNAMatched = -1 
    bestShift = -1
    bestNumCompared = -1
    
    for i in range(len(globalReference)):
        shift,counter,numCompared = getAlignmentsForBestMatch(globalReferenceScaffold[i],myRead)[0]
        
        if(counter >bestNumDNAMatched):
            bestBG = i
            bestNumDNAMatched = counter
            bestNumCompared=numCompared
            bestShift=shift
            
    return bestBG,(bestShift,bestNumDNAMatched,bestNumCompared)
        

def getDNAAndAA(forwardRecord,reverseRecord):
    
    ##################################
    alignToRef = True;
    numMutations = globalNumMutations
    numMaxOverlap = globalNumMaxOverlap
    ##################################
    
    forwardReads = str(forwardRecord.seq)
    reverseReads = str(reverseRecord.seq.reverse_complement())
    
    referenceShift,forwardShift,reverseShift="","",""
    
    refIndex,alignForwardToRef = getBestBG(forwardReads) 
    # To get the best background, all references have to be tested. 
    # Pull out the alignment to use in the rest of the process
    
    alignReverseToRef=None
    reference = globalReferenceScaffold[refIndex]
    
    if(alignForwardToRef[0]<0): #Alignment to reference should be moved up
        referenceShift = "."*-alignForwardToRef[0]
        forwardShift = ""
    else:
        referenceShift=""
        forwardShift=" "*alignForwardToRef[0]


    ##Typically, match the reverse read to the tail end of forward read
    alignReverseToForward = getAlignmentsForTailMatch(forwardReads,
                                                      reverseReads,
                                                      maxOverlap=numMaxOverlap,
                                                      allowableMutations=numMutations,
                                                      verbose=False)[0]
    reverseShift=forwardShift+" "*alignReverseToForward[0]
    
    #MAP BOTH TO REFERENCE If there is <1 overlap
    if (alignReverseToForward[0]>=len(forwardReads)-globalMinNumOverlap):
        alignReverseToRef = getAlignmentsForBestMatch(reference,reverseReads)[0]
        reverseShift=referenceShift+" "*alignReverseToRef[0]
        
    string0 = forwardShift + forwardReads
    string1 = reverseShift+ reverseReads

    #Assemble the transcript
    
    forward_Quality = forwardRecord.letter_annotations['phred_quality']
    reverse_Quality = reverseRecord.letter_annotations['phred_quality'][::-1]
    
    alignedFQuality = [0]*len(forwardShift)
    alignedFQuality.extend(forward_Quality)
    alignedRQuality = [0]*len(reverseShift)
    alignedRQuality.extend(reverse_Quality)

    
    myRead = []
    myReadQuality=[]
    for i in range(len(string1)):
        
        # Append nt to string if there is no reverse read
        if(i<len(string0)):
            if(string0[i]=='.'):
                pass
            else:
                if(string1[i]==' '):
                    myRead.append(string0[i])
                    myReadQuality.append(alignedFQuality[i])
                else:
                    if(alignedFQuality[i]>=alignedRQuality[i]):
                        myRead.append(string0[i])
                        myReadQuality.append(alignedFQuality[i])

                    elif(alignedFQuality[i]<alignedRQuality[i]):
                        myRead.append(string1[i])
                        myReadQuality.append(alignedRQuality[i])
                    else:
                        myRead.append('N')
                        myReadQuality.append(0)

        else:
            if(string1[i]==' '):
                myRead.append('N')
                myReadQuality.append(0)

            else:
                myRead.append(string1[i])
                myReadQuality.append(alignedRQuality[i])

    

    
    
    
    myDNA = "".join(myRead)
    
    myFrameShift = len(referenceShift)
    myExtractedDNA = myDNA[(globalRangeToExtract[0]+myFrameShift):(globalRangeToExtract[1]+myFrameShift)]
    myExtractedDNA = myExtractedDNA.replace(" ","N")
    myExtractedAA =  Seq.Seq(myExtractedDNA).translate()

    myExtractedQuality =myReadQuality[(globalRangeToExtract[0]+myFrameShift):(globalRangeToExtract[1]+myFrameShift)]
    #if(refIndex==1):
    
    #For debugging purposes - Uhash these 3 lines
    #print( 'R',referenceShift+reference)                  
    #print('+',forwardShift + forwardReads)
    #print('-',reverseShift+ reverseReads)
    
    #print(myExtractedAA)
    #print(myExtractedQuality)
    #print('---------------------------------')
    #print(alignedFQuality)
    #print(alignedRQuality)
        #print str("".join(myRead)).rjust(120)

    print (myExtractedAA)
    return refIndex,myExtractedDNA, myExtractedAA, myExtractedQuality

if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
#processFile("/home/vxue/luther/SORTCERY/SORTCERY2_all/barcode_",'x',"ref","expected")
