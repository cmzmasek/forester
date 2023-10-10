CD_HIT="/Users/czmasek/anaconda3/bin/cd-hit"

if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <cutoff> <indir> <outdir>" >&2
  exit 1
fi

cutoff=$1
indir=$2
outdir=$3

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
            echo "Executing: $CD_HIT -c $cutoff -i $indir/$name.fasta -o $outdir/$name.fasta:"
            $CD_HIT -c $cutoff -i $indir/$name.fasta -o $outdir/$name.fasta
        fi
    fi
done