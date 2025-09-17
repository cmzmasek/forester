re="(.+)_ni\.fasta$"
#re="(.+)\.fasta$"
MAFFT="/Users/czmasek/anaconda3/bin/mafft"

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <dir>" >&2
  exit 1
fi

indir=$1

if [ -z "$(ls -A $indir)" ]; then
  echo "Directory is empty"
  exit 0
fi

cutoff=3
for i in $indir/*; do
  if test -f "$i"; then
    if [[ $i =~ $re ]]; then
      name=${BASH_REMATCH[1]}
      echo $name
      count=$(grep -c '>' $i)
      echo $count
      if [ $count -gt $cutoff ]; then
        $MAFFT --maxiterate 1000 --localpair --thread 8 $i >${name}_mafft_1000_l.fasta
        rc=$?
        if [[ $rc != 0 ]]; then
          exit $rc
        fi
      fi
    fi
  fi
done
