#!/bin/bash
cd seqframe;
mkdir summary;
for each in `ls seqframe_*`;
    do sort $each | uniq -c > summary/$each.smry;
done
