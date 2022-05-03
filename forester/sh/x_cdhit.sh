CD_HIT="/home/zma/SOFTWARE/CDHIT/cdhit-4.6.6/cd-hit"
CUTOFF=0.98

indir=$1
outdir=$2

if [ ! -d "$outdir" ]; then
    mkdir $outdir
fi

if [ ! -z "$(ls -A $outdir)" ]; then
    echo "Output directory is not empty"
    exit 0
fi

if [ -z "$(ls -A $indir)" ]; then
    echo "Input directory is empty"
    exit 0
fi

fasta_re="/(.+)\.fasta$"

for i in $indir/* 
do
    if test -f "$i" 
    then
        if [[ $i =~ $fasta_re ]]
        then
            name=${BASH_REMATCH[1]}
            echo "Executing: $CD_HIT -c $CUTOFF -i $indir/$name.fasta -o $outdir/$name.fasta:"
            $CD_HIT -c $CUTOFF -i $indir/$name.fasta -o $outdir/$name.fasta
        fi
    fi
done