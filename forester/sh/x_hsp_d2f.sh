re="(.+)_hmmscan$"
for i in * 
do
    if test -f "$i" 
    then
        if [[ $i =~ $re ]]
        then
            name=${BASH_REMATCH[1]}
            echo $name
            ruby /home/zma/git/forester/forester/ruby/evoruby/exe/hsp.rb $i
            rc=$?
            if [[ $rc != 0 ]]
            then
                exit $rc
            fi
            ruby /home/zma/git/forester/forester/ruby/evoruby/exe/d2f.rb -o ${name}_hmmscan_domain_table
            rc=$?
            if [[ $rc != 0 ]]
            then
                exit $rc
            fi
        fi
    fi
done
