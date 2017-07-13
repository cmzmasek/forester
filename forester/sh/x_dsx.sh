re="(.+)_hmmscan$"
for i in * 
do
    if test -f "$i" 
    then
        if [[ $i =~ $re ]]
        then
            name=${BASH_REMATCH[1]}
            echo $name
            ruby /home/zma/git/forester/forester/ruby/evoruby/exe/dsx.rb -d -e=1e-6 -l=50 ${name} ${name}_hmmscan
            rc=$?
            if [[ $rc != 0 ]]
            then
                exit $rc
            fi
        fi
    fi
done