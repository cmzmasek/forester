CD_HIT="/usr/bin/cd-hit"
CUTOFF=0.9995
MAX=1000
FASTA_RE="/(.+)\.fasta$"

indir=$1
outdir=$2

if [ ! -d "$outdir" ]; then
    mkdir $outdir
fi

if [ -z "$(ls -A $indir)" ]; then
    echo "Input directory is empty"
    exit 0
fi

for i in $indir/* 
do
    if test -f "$i" 
    then
        if [[ $i =~ $FASTA_RE ]]
        then
            name=${BASH_REMATCH[1]}
            echo $name
            if [ -s $outdir/$name.fasta ]
            then
               echo $name already done
            else
                count=$(grep -c '>' $i)
                echo $count
                if [ $count -gt $MAX ]
                then
                    echo "Executing: $CD_HIT -c $CUTOFF -T 0 -M 2000 -i $indir/$name.fasta -o $outdir/$name.fasta:"
                    $CD_HIT -c $CUTOFF -T 0 -M 2000 -i $indir/$name.fasta -o $outdir/$name.fasta
                else
                    echo "Executing: cp $indir/$name.fasta $outdir/:"
                    cp $indir/$name.fasta $outdir/    
                fi
            fi       
        fi
    fi
done
