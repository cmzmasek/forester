re="(.+)_mafft_1000_g_05_20"
for i in * 
do
    if test -f "$i" 
    then
        if [[ $i =~ $re ]]
        then
            name=${BASH_REMATCH[1]}
            echo $name
            perl /home/zma/git/forester/forester/perl/phylo_pl.pl -B1000Wq@1S9X $i ${name}_mafft_1000_g_05_20_tree
            rc=$?
            if [[ $rc != 0 ]]
            then
                exit $rc
            fi
        fi
    fi
done