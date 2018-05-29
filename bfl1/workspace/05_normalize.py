import pandas as pd
import numpy as np

dirName = "/home/vxue/data/SORTCERY_PUBLICATION/bfl1/workspace/"
experimentNames = "/home/vxue/data/SORTCERY_PUBLICATION/bfl1/experimentNames.txt"



normalization = pd.read_csv(experimentNames,header=None,names=['name']+list(map(str,range(12))))


def getNormlization(row,columnSums,normRate):
    reads = row[list(range(12))]
    #print(reads)
    percentage = reads / columnSums 
    
    cells = np.multiply(percentage,normRate)
    
    return  cells / np.sum(cells)

for suffix,fileOutput in zip(["","_dna"],['seqframe','dnaframe']):
    for idx,row in normalization.iterrows():
        
        i = row['name']
        rate = row[list(map(str,list(range(12))))]
        
        scaffold0SeqsFile = dirName+i+suffix+"/unique_seqs_"+fileOutput+"_0"
        scaffold1SeqsFile = dirName+i+suffix+"/unique_seqs_"+fileOutput+"_1"
        
        
        scaffold0CountFile = dirName+i+suffix+"/cell_counts_"+fileOutput+"_0"
        scaffold1CountFile = dirName+i+suffix+"/cell_counts_"+fileOutput+"_1"
        
        
        #print(scaffold1SeqsFile)
        #print(scaffold1CountFile)
        
        
        df0 = pd.read_csv(scaffold0CountFile,delimiter="\s+",header=None)
        df0['seq'] = pd.read_csv(scaffold0SeqsFile,header=None)[0]
        df0['bg'] = 0
        
        df1 = pd.read_csv(scaffold1CountFile,delimiter="\s+",header=None)
        df1['seq'] = pd.read_csv(scaffold1SeqsFile,header=None)[0]
        df1['bg'] = 1
        
        concat = pd.concat([df0,df1])
        concat['tot'] = concat[list(range(12))].sum(axis=1)
        
        
        subset = concat[concat.tot>=100]
        
        sums = subset[list(range(12))].sum(axis=0)
        
        normed = subset.apply(lambda x: getNormlization(x,sums,rate),axis=1)
        
        subset[['bg','seq']].to_csv(dirName+i+suffix+"/sequences_cutoff_100",header=None,index=None,sep=' ')
        normed.to_csv(dirName+i+suffix+"/norm_distributions_cutoff_100",header=None,index=None,sep=' ',float_format='%.6f')
        
        print(i,suffix,len(normed))
