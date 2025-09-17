VERSION="1.0.1"
TAP="/usr/bin/ruby /Users/czmasek/IdeaProjects/forester/forester/ruby/evoruby/exe/tap.rb"
MSA_PRO="/usr/bin/ruby /Users/czmasek/IdeaProjects/forester/forester/ruby/evoruby/exe/msa_pro.rb"
DECORATOR="/usr/bin/ruby /Users/czmasek/IdeaProjects/forester/forester/ruby/evoruby/exe/phylogenies_decorator.rb"
PHYLO_PL="perl /Users/czmasek/IdeaProjects/forester/forester/perl/phylo_pl.pl"
MAFFT="/Users/czmasek/anaconda3/bin/mafft"
MSA_RENAME="python3 /Users/czmasek/Dropbox/PROG/PYTHON/PYCHARM_PROJECTS/TWO/msa_rename.py"

if [ "$#" -ne 5 ]; then
  echo "Usage: mktree.sh <mafft options> <msa_pro options> <phylopl options> <input suffix> <workdir>" >&2
  echo "    Example 1: mktree.sh \"--auto --thread -1\" \"\" \"-B1Zq\" .fasta trees"
  echo "    Example 2: mktree.sh \"--retree 1 --thread -1\" \"\" \"-B1Wnq\" .fasta trees"
  echo "    Example 3: mktree.sh \"--maxiterate 1000 --localpair --thread 8\" \"-rr=0.5 -rsgr=0.5\" \"-B100Znq\" .fasta trees"
  echo "    Example 4: mktree.sh \"--maxiterate 1000 --localpair --thread 8\" \"-rr=0.5 -rsgr=0.5\" \"-B100Wnqxg\" .fasta trees"
  exit 1
fi

mafft_opts=$1
msapro_opts=$2
phylopl_opts=$3
input_suffix=$4
workdir=$5

echo "Version         : $VERSION"
echo "MAFFT options   : $mafft_opts"
echo "msa_pro options : $msapro_opts"
echo "phylo.pl options: $phylopl_opts"
echo "Input suffix    : $input_suffix"
echo "Work dir        : $workdir"
echo ""
echo ""

if [ ! -d "$workdir" ]; then
  echo "Working directory does not exist"
  exit 0
fi

if [ -z "$(ls -A $workdir)" ]; then
  echo "Working directory is empty"
  exit 0
fi

input_suffix_re=".*/(.+)$input_suffix$"

for i in $workdir/*; do
  if test -f "$i"; then
    if [[ ! $i =~ "_ni.fasta" && ! $i =~ "_mafft.fasta" && $i =~ $input_suffix_re ]]; then
      name=${BASH_REMATCH[1]}

      count=$(grep -c '>' $i)
      echo "    Count for $name$input_suffix $count"
      if [ $count -gt 3 ]; then

        echo "    Working on: $name$input_suffix"

        if [ ! -f $workdir/${name}_ni.fasta ]; then
          echo "        Executing: $TAP $workdir/$name$input_suffix:"
          $TAP $workdir/$name$input_suffix
          rc=$?
          if [[ $rc != 0 ]]; then
            exit $rc
          fi
        fi

        if [ ! -f $workdir/${name}_ni_mafft.fasta ]; then
          echo "        Executing: $MAFFT $mafft_opts $workdir/${name}_ni.fasta > $workdir/${name}_ni_mafft.fasta:"
          $MAFFT $mafft_opts $workdir/${name}_ni.fasta >$workdir/${name}_ni_mafft.fasta
          rc=$?
          if [[ $rc != 0 ]]; then
            exit $rc
          fi
        fi

        if [ ! -f $workdir/${name}_mafft.fasta ]; then
          echo "        Executing: $MSA_RENAME $workdir/${name}_ni_mafft.fasta $workdir/${name}.nim $workdir/${name}_mafft.fasta:"
          $MSA_RENAME $workdir/${name}_ni_mafft.fasta $workdir/${name}.nim $workdir/${name}_mafft.fasta
          rc=$?
          if [[ $rc != 0 ]]; then
            exit $rc
          fi
        fi

        if [ ! -f $workdir/${name}_ni_mafft ]; then
          echo "        Executing: $MSA_PRO -i=f -o=p -d -c $msapro_opts  $workdir/${name}_ni_mafft.fasta $workdir/${name}_ni_mafft:"
          $MSA_PRO -i=f -o=p -d -c $msapro_opts $workdir/${name}_ni_mafft.fasta $workdir/${name}_ni_mafft
          rc=$?
          if [[ $rc != 0 ]]; then
            exit $rc
          fi
        fi

        if [ ! -f $workdir/${name}_ni_mafft_tree.log ]; then
          echo "        Executing: $PHYLO_PL -$phylopl_opts $workdir/${name}_ni_mafft $workdir/${name}_ni_mafft_tree:"
          $PHYLO_PL -$phylopl_opts $workdir/${name}_ni_mafft $workdir/${name}_ni_mafft_tree
          rc=$?
          if [[ $rc != 0 ]]; then
            exit $rc
          fi
        fi
      fi
    fi
  fi
done

if [ ! -d "$workdir/deco" ]; then
  mkdir $workdir/deco
fi
rm $workdir/deco/*
cp $workdir/*.xml $workdir/deco/
cp $workdir/*.nim $workdir/deco/
cp $workdir/*_ni.fasta $workdir/deco/

current_dir=$(pwd)
cd $workdir/deco/
echo "        Executing: $DECORATOR -nd -v .xml _d.xml:"
$DECORATOR -nd -v .xml _d.xml
rc=$?
if [[ $rc != 0 ]]; then
  exit $rc
fi
cd $current_dir

rm $workdir/deco/*.nim
rm $workdir/deco/*_ni.fasta

rm $workdir/*_fme_stats.txt

echo "Version         : $VERSION"
echo "MAFFT options   : $mafft_opts"
echo "msa_pro options : $msapro_opts"
echo "phylo.pl options: $phylopl_opts"
echo "Input suffix    : $input_suffix"
echo "Work dir        : $workdir"
echo "mktree successfully completed"
echo ""
