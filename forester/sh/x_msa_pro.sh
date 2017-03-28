re="(.+)_mafft_1000_g\.aln"
for i in * 
do
    if test -f "$i" 
    then
        if [[ $i =~ $re ]]
        then
            name=${BASH_REMATCH[1]}
            echo $name
            ruby /home/zma/git/forester/forester/ruby/evoruby/exe/msa_pro.rb -i=f -o=p -d -c -rr=0.5 -rsl=20 $i ${name}_mafft_1000_g_05_20
            rc=$?
            if [[ $rc != 0 ]]
            then
                exit $rc
            fi
        fi
    fi
done