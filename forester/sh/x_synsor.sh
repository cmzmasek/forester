SYNSOR=""

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <indir>" >&2
  exit 1
fi

indir=$1

if [ -z "$(ls -A $indir)" ]; then
    echo "Directory is empty"
    exit 0
fi

re="(.+)\.fasta$"

for i in $indir/*
do
    if test -f "$i"
    then
        if [[ $i =~ $re ]]
        then
            name=${BASH_REMATCH[1]}
            echo $name
            echo SYNSOR  -i $name.fasta -k 7 -m model -o ${name}_results
            rc=$?
            if [[ $rc != 0 ]]
            then
                exit $rc
            fi
        fi
    fi
done