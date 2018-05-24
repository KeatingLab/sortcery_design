#!/home/vxue/anaconda/bin/python

#Tally up all the observed nucleotides at each site.

import numpy as np
import pickle 
dnaIndex = 'ATCGN'
myTally=[None,None]
for bg in range(2):
    
    myTally[bg] = np.zeros((57,5))
    for barcode in range(168):
        print bg,barcode
        myFile = '/home/vxue/luther/SORTCERY/SORTCERY_justin/dnaframe/dnaframe_'+str(bg)+'_code_'+str(barcode)
        with open(myFile,'r') as readStream:
            for each in readStream:
                for i in range(57):
                    myTally[bg][i][dnaIndex.index(each[i])]+=1
                    
pickle.dump(myTally,open("/home/vxue/luther/SORTCERY/SORTCERY_justin/dnaCounts.pickle",'w'))
