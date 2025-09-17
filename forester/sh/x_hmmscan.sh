HMMSCAN="/usr/local/bin/hmmscan"
PFAM="/Users/czmasek/DATA/PFAM_380X/Pfam-A.hmm"
RE="(.+)_ni\.fasta$"

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <indir>" >&2
  exit 1
fi

indir=$1

if [ -z "$(ls -A $indir)" ]; then
  echo "Directory is empty"
  exit 0
fi

for i in $indir/*; do
  if test -f "$i"; then
    if [[ $i =~ $RE ]]; then
      name=${BASH_REMATCH[1]}
      echo $name
      $HMMSCAN --domtblout ${name}_hmmscan -E 20 --domE 20 --noali $PFAM ${name}_ni.fasta
      rc=$?
      if [[ $rc != 0 ]]; then
        exit $rc
      fi
    fi
  fi
done
