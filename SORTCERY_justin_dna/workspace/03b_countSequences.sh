#!/bin/bash
cd dnaframe;
mkdir summary;
for each in `ls dnaframe_*`;
    do sort $each | uniq -c > summary/$each.smry;
done
