re="(.+)-Alignment\.txt$"
for i in *
do
    if test -f "$i" 
    then
        if [[ $i =~ $re ]]
        then
            name=${BASH_REMATCH[1]}
            echo $name
            java -Xmx8048m -cp /home/lambda/git/forester/forester/java/forester.jar org.forester.application.blast_pars $i > ${name}.fasta
            rc=$?
            if [[ $rc != 0 ]]
            then
                exit $rc
            fi
        fi
    fi
done