# x_phylo
# -------

# Version: 1.00
# Last modified: 180504

# Usage: x_phylo.sh <indir> <outdir>


# OPTIONS
# -------

MAFFT_OPTIONS="--maxiterate 1000 --localpair"
#MAFFT_OPTIONS="--auto"
MSA_PRO_OPTIONS="-rr=0.5 -rsl=20"
PHYLO_PL_OPTIONS="-B100Wq@1S9X"


# PROGRAMS
# --------

TAP="/home/zma/git/forester/forester/ruby/evoruby/exe/tap.rb"
MSA_PRO="/home/zma/git/forester/forester/ruby/evoruby/exe/msa_pro.rb"
MAFFT="/usr/local/bin/mafft"
PHYLO_PL="/home/zma/git/forester/forester/perl/phylo_pl.pl"


# PRE
# ---

if [[ $# -ne 2 ]] ; then
    echo "Usage: x_phylo.sh <indir> <outdir>"
    exit 0
fi

indir=$1
outdir=$2

echo "Input directory  : " $indir
echo "Output directory : " $outdir
echo "MAFFT options    : " $MAFFT_OPTIONS
echo "MSA PRO options  : " $MSA_PRO_OPTIONS
echo "PHYLO PL options : " $PHYLO_PL_OPTIONS

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

mkdir $outdir/maps
mkdir $outdir/seqs
mkdir $outdir/msas
mkdir $outdir/phylo_trees
mkdir $outdir/phylo_mlt
mkdir $outdir/phylo_aux
mkdir $outdir/phylo_logs


# TAP
# ---

fasta_re="(.+)\.fasta$"

for i in $indir/* 
do
    if test -f "$i" 
    then
        if [[ $i =~ $fasta_re ]]
        then
            name=${BASH_REMATCH[1]}
            echo "Working on " $name
            ruby $TAP -t $i
            rc=$?
            if [[ $rc != 0 ]]
            then
                exit $rc
            fi
        fi
    fi
done

mv $indir/*.nim $outdir/maps
mv $indir/*_ni.fasta $outdir/seqs


# MAFFT
# -----

ni_fasta_re="(.+)_ni\.fasta$"
cutoff=3

for i in $outdir/seqs/* 
do
    if test -f "$i" 
    then
        if [[ $i =~ $ni_fasta_re ]]
        then
            name=${BASH_REMATCH[1]}
            echo "Working on: " $name
            count=$(grep -c '>' $i)
            echo $count
            if [ $count -gt $cutoff ]
            then
                $MAFFT $MAFFT_OPTIONS $i > ${name}_mafft.fasta
                rc=$?
                if [[ $rc != 0 ]]
                then
                    exit $rc
                fi
            fi
        fi
    fi
done

mv $outdir/seqs/*_mafft.fasta $outdir/msas/


# MSA-PRO
# ------

mafft_out_re="(.+)_mafft\.fasta$"
for i in $outdir/msas/*  
do
    if test -f "$i" 
    then
        if [[ $i =~ $mafft_out_re ]]
        then
            name=${BASH_REMATCH[1]}
            echo $name
            ruby $MSA_PRO -i=f -o=p -d -c $MSA_PRO_OPTIONS $i ${name}_mafft_p
            rc=$?
            if [[ $rc != 0 ]]
            then
                exit $rc
            fi
        fi
    fi
done


# PHYLO_PL
# --------

msa_pro_out_re="(.+)_mafft_p$"
for i in $outdir/msas/* 
do
    if test -f "$i" 
    then
        if [[ $i =~ $msa_pro_out_re ]]
        then
            name=${BASH_REMATCH[1]}
            echo $name
            perl $PHYLO_PL $PHYLO_PL_OPTIONS $i ${name}_mafft_tree
            rc=$?
            if [[ $rc != 0 ]]
            then
                exit $rc
            fi
            rm ${name}_mafft_tree.mpwd 
        fi
    fi
done

mv $outdir/msas/*_tree_*.xml $outdir/phylo_trees/
mv $outdir/msas/*_tree.log $outdir/phylo_logs/
mv $outdir/msas/*_tree_puzzle_outfile $outdir/phylo_aux/
mv $outdir/msas/*_tree_*.mlt $outdir/phylo_mlt/


# cd $outdir/phylo_trees/
# phylogenies_decorator -nd -ns  -tc .xml _d.xml ../maps/
# cd ..
# mkdir phylo_trees_decorated
# mv phylo_trees/00_phylogenies_decorator.log phylo_trees_decorated
# mv phylo_trees/*_d.xml phylo_trees_decorated
# gsdi -g -r -t -s=.xml phylo_trees_decorated/ ../../Coronaviridae_Taxonomy2.xml phylo_trees_gsdi
# gsdi -g -r -R -t -s=.xml phylo_trees_decorated/ ../../Coronaviridae_Taxonomy2.xml phylo_trees_gsdi_R

