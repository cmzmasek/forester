re="Q_(.+)\.fasta"

for i in * 
do
    if test -f "$i" 
    then
        if [[ $i =~ $re ]]
        then
            name=${BASH_REMATCH[1]}
            echo $name
            /usr/local/bin/mafft --keeplength --add $i CDS_aligned_partials_MAFFT_merged_ni.fasta > REF_PLUS_Q_${name}.fasta
            rc=$?
            if [[ $rc != 0 ]]
            then
                exit $rc
            fi
            
           /home/zma/SOFTWARE/PHYLO/PPLACER/PPLACER11ALPHA19/pplacer-Linux-v1.1.alpha19/pplacer -m GTR -t RAxML_bestTree.CDS_aligned_partials_MAFFT_merged_GTR -s RAxML_info.CDS_aligned_partials_MAFFT_merged_GTR REF_PLUS_Q_${name}.fasta -o PP_Q_${name}.json
           rc=$?
           if [[ $rc != 0 ]]
           then
               exit $rc
           fi
            
           /home/zma/SOFTWARE/PHYLO/PPLACER/PPLACER11ALPHA19/pplacer-Linux-v1.1.alpha19/guppy sing PP_Q_${name}.json
           rc=$?
           if [[ $rc != 0 ]]
           then
               exit $rc
           fi
        fi
    fi
done