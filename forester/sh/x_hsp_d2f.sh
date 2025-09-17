HSP="/Users/czmasek/IdeaProjects/forester/forester/ruby/evoruby/exe/hsp.rb"
D2F="/Users/czmasek/IdeaProjects/forester/forester/ruby/evoruby/exe/d2f.rb"
RE="(.+)_hmmscan$"

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
      ruby $HSP $i
      rc=$?
      if [[ $rc != 0 ]]; then
        exit $rc
      fi
      ruby $D2F -o ${name}_hmmscan_domain_table ${name}_ni.fasta
      rc=$?
      if [[ $rc != 0 ]]; then
        exit $rc
      fi
    fi
  fi
done
