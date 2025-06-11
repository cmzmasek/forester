VERSION="1.0.0"
CLEAN_FASTA="python /Users/czmasek/Dropbox/PROG/PYTHON/PYCHARM_PROJECTS/TWO/clean_fasta.py"
CDHIT="/Users/czmasek/anaconda3/bin/cd-hit"
MK_TREE="bash /Users/czmasek/IdeaProjects/forester/forester/sh/mktree.sh"
VIPR_X="java -Xmx8048m -cp /Users/czmasek/IdeaProjects/forester/forester/java/forester.jar org.forester.application.vipr_x6"


genomes_file=$1
annotation_file=$2
seq_type=$3
minimal_length_a=$4
minimal_length_b=$5
char_ratio=$6
clustering_cutoff_a=$7
clustering_cutoff_b=$8
mafft_opts=$9
tree_inference_opt=$10
date_str=$
tree_a_desc=$
tree_b_desc=$
tree_a_output_name=
tree_b_output_name=
workdir=$

if [ "$#" -ne  ]; then
  echo "Usage: " >&2
  echo "    Example 1: ref_phylo_pl.sh ref_phylo_pl.sh BVBRC_1900_2025.fasta ../JUNE_10_2025/BVBRC_genome.txt na 12000 450 0.999 0.95 workdir"
   exit 1
fi

echo "Version             : $VERSION"
echo "Genomes file        : $genomes_file"
echo "Annotations file    : $annotation_file"
echo "Sequence type       : $seq_type"
echo "Min length A        : $minimal_length_a"
echo "Min length B        : $minimal_length_b"
echo "Ratio               : $char_ratio"
echo "Clustering cutoff A : $clustering_cutoff_a"
echo "Clustering cutoff B : $clustering_cutoff_b"
echo "Date string         : $date_str"
echo "Tree desc A         : $tree_a_desc"
echo "Tree desc B         : $tree_b_desc"
echo "Work dir            : $workdir"
echo ""
echo ""

if [ ! -d "$workdir" ]; then
  mkdir $workdir
fi

if [ ! -z "$(ls -A $workdir)" ]; then
  echo "Working directory is not empty"
  exit -1
fi

if [ ! -z "$(ls -A $workdir)" ]; then
  echo "Tree already exists"
  exit -1
fi

if [ ! -z "$(ls -A $workdir)" ]; then
  echo "Tree already exists"
  exit -1
fi

$CLEAN_FASTA -ml $minimal_length_a -r $char_ratio -t $seq_type $genomes_file $workdir/genomes_clean_a.fasta

$CLEAN_FASTA -ml $minimal_length_b -r $char_ratio -t $seq_type $genomes_file $workdir/genomes_clean_b.fasta

$CDHIT -c $clustering_cutoff_a -d 0 -g 1 -i $workdir/genomes_clean_a.fasta -o $workdir/genomes_clean_a_cdhit.fasta

$CDHIT -c $clustering_cutoff_b -d 0 -g 1 -i $workdir/genomes_clean_b.fasta -o $workdir/genomes_clean_b_cdhit.fasta

mkdir $workdir/trees

cp $workdir/genomes_clean_a_cdhit.fasta $workdir/trees/
cp $workdir/genomes_clean_b_cdhit.fasta $workdir/trees/

$MK_TREE "--auto --thread 8" "-rr=0.2 -rsl=450" "-B100Zxm" .fasta $workdir/trees/

$VIPR_X $workdir/trees/deco/genomes_clean_a_cdhit_ni_mafft_tree_raxml_d.xml $annotation_file measles_virus_compl_genome_cdhit099.xml "M. hominis genome 2025-XX-XX cdhit0.99 RAxML GTR G4"

