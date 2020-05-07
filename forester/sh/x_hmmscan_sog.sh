#!/bin/bash
re="(.+)\.fasta$"
for i in * 
do
    if test -f "$i" 
    then
        if [[ $i =~ $re ]]
        then
            name=${BASH_REMATCH[1]}
            echo $name
            /home/lambda/SOFT/hmmer-3.3/src/hmmscan --nobias --domtblout ${name}.hmmscan -E 10 --domE 10 --noali ~/DATA/PFAM33/Pfam-A.hmm ${name}.fasta
        fi
    fi
done

