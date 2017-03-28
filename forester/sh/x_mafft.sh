re="(.+)_ni\.fasta"
cutoff=3
for i in * 
do
    if test -f "$i" 
    then
        if [[ $i =~ $re ]]
        then
            name=${BASH_REMATCH[1]}
            echo $name
            count=$(grep -c '>' $i)
            echo $count
            if [ $count -gt $cutoff ]
            then
                /usr/local/bin/mafft --maxiterate 1000 --globalpair $i > ${name}_mafft_1000_g.aln
                rc=$?
                if [[ $rc != 0 ]]
                then
                    exit $rc
                fi
            fi
        fi
    fi
done