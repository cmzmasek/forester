#!/usr/bin/env bash
re="(.+__.+)\$"
for i in *
do
    if test -f "$i"
    then
        if [[ $i =~ $re ]]
        then
            name=${BASH_REMATCH[1]}
            echo $name
            /usr/local/bin/hmmscan --domtblout ${name}.hmmscan -E 1 --domE 1 --noali ~/DATA/PFAM340X/Pfam-A.hmm ${name} >/dev/null &
        fi
    fi
done
