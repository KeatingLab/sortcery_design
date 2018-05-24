#!/bin/bash
for each in `sed -n "1~13 p" ../multiplex.txt`;
    do 
    myDir=$each"_dna";
    echo $myDir;
    /home/vxue/scratch/Applications/usearch8.0.1517_i86linux32 -cluster_otus $myDir/diversity.fasta -sizein -sizeout -otu_radius_pct 3 -otus $myDir/rep.fasta -uparseout $myDir/out.up
done;

