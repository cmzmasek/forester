re="(.+)_ni\.fasta$"
#re="(.+)\.fasta$"
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
                /usr/local/bin/mafft --maxiterate 1000 --localpair $i > ${name}_mafft_1000_l.fasta
                rc=$?
                if [[ $rc != 0 ]]
                then
                    exit $rc
                fi
            fi
        fi
    fi
done