re="(.+)_mafft_1000_l\.fasta$"

MSA_PRO="/Users/czmasek/IdeaProjects/forester/forester/ruby/evoruby/exe/msa_pro.rb"

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <dir>" >&2
  exit 1
fi

indir=$1

if [ -z "$(ls -A $indir)" ]; then
  echo "Directory is empty"
  exit 0
fi

for i in $indir/*; do
  if test -f "$i"; then
    if [[ $i =~ $re ]]; then
      name=${BASH_REMATCH[1]}
      echo $name
      ruby $MSA_PRO -i=f -o=p -d -c $i ${name}_mafft_1000_l
      rc=$?
      if [[ $rc != 0 ]]; then
        exit $rc
      fi
    fi
  fi
done
