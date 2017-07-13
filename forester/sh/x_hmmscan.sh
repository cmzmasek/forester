re="(.+)_ni\.fasta$"
for i in * 
do
    if test -f "$i" 
    then
        if [[ $i =~ $re ]]
        then
            name=${BASH_REMATCH[1]}
            echo $name
            /usr/local/bin/hmmscan --max --domtblout ${name}_hmmscan -E 20 --domE 20 --noali ~/DATA/PFAM/PFAM_30/Pfam-A.hmm ${name}_ni.fasta
        fi
    fi
done