re="(.+)_mafft_1000_l$"

PHYLO_PL="/Users/czmasek/IdeaProjects/forester/forester/perl/phylo_pl.pl"

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <dir>" >&2
  exit 1
fi

indir=$1

if [ -z "$(ls -A $indir)" ]; then
    echo "Directory is empty"
    exit 0
fi

for i in $indir/*
do
    if test -f "$i" 
    then
        if [[ $i =~ $re ]]
        then
            name=${BASH_REMATCH[1]}
            echo $name
            perl $PHYLO_PL -B100Wq $i ${name}_mafft_1000_l_tree_WAG
            rc=$?
            if [[ $rc != 0 ]]
            then
                exit $rc
            fi
        fi
    fi
done