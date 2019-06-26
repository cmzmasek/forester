re="(.+)\.fasta$"
for i in * 
do
    if test -f "$i" 
    then
        if [[ $i =~ $re ]]
        then
            name=${BASH_REMATCH[1]}
            echo $name
            ruby /home/zma/git/forester/forester/ruby/evoruby/exe/tap.rb -t $i
            rc=$?
            if [[ $rc != 0 ]]
            then
                exit $rc
            fi
        fi
    fi
done

I had high hopes for this movie. That is why I watched it on Thursday night. Unfortunately, it turned out a complete disappointment. The story is boring and predictable, the action sequences amateurish in their execution, and the humor falls flat. The worst of all is the main actress, though. What were they thinking by casting her? She does not fit the role at all! She has zero charisma and actually seemed bored during filming.