#!/bin/bash
for file_name in $1/*; do
    echo $file_name
    /home/czmasek/SOFTWARE/GENSCAN/genscan /home/czmasek/SOFTWARE/GENSCAN/HumanIso.smat $file_name -cds >> $2
done

