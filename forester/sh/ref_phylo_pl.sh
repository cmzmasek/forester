VERSION="0.0.1"
CLEAN_FASTA="python /Users/czmasek/Dropbox/PROG/PYTHON/PYCHARM_PROJECTS/TWO/clean_fasta.py"
CDHIT="/Users/czmasek/anaconda3/bin/cd-hit"
MK_TREE="bash /Users/czmasek/IdeaProjects/forester/forester/sh/mktree.sh"
VIPR_X="java -Xmx8048m -cp /Users/czmasek/IdeaProjects/forester/forester/java/forester.jar org.forester.application.vipr_x6"

if [ "$#" -ne 8 ]; then
  echo "Usage: " >&2
  echo "    Example 1: ref_phylo_pl.sh ref_phylo_pl.sh BVBRC_1900_2025.fasta ../JUNE_10_2025/BVBRC_genome.txt na 12000 450 0.999 0.95 workdir"
  exit 1
fi

options_file=$1
genomes_file=$2
annotation_file=$3
tree_a_desc=$4
tree_b_desc=$5
tree_a_output_name=$6
tree_b_output_name=$7
workdir=$8

if [ ! -f "$options_file" ]; then
  echo "$options_file does not exist"
  exit -1
fi

seq_type="?"
minimal_length_a=-1
char_ratio=-1
for line in $(cat $options_file); do
  echo $line
  if [[ "$line" =~ SEQ_TYPE=(.+) ]]; then
    seq_type="${BASH_REMATCH[1]}"
  fi
  if [[ "$line" =~ MINIMAL_LENGTH_A=(.+) ]]; then
    minimal_length_a="${BASH_REMATCH[1]}"
  fi
  if [[ "$line" =~ MINIMAL_LENGTH_B=(.+) ]]; then
    minimal_length_b="${BASH_REMATCH[1]}"
  fi
  if [[ "$line" =~ CHAR_RATIO=(.+) ]]; then
    char_ratio="${BASH_REMATCH[1]}"
  fi

done

echo "Version             : $VERSION"
echo "Options file        : $options_file"
echo "Genomes file A      : $genomes_file_a"
echo "Genomes file B      : $genomes_file_b"
echo "Annotations file    : $annotation_file"
echo "Sequence type       : $seq_type"
echo "Min length A        : $minimal_length_a"
echo "Min length B        : $minimal_length_b"
echo "Ratio               : $char_ratio"
echo "Clustering cutoff A : $clustering_cutoff_a"
echo "Clustering cutoff B : $clustering_cutoff_b"
echo "MAFFT options A     : $mafft_opts_a"
echo "MAFFT options B     : $mafft_opts_b"
echo "MSA pro options A   : $msa_pro_opts_a"
echo "MSA pro options B   : $msa_pro_opts_b"
echo "Tree inf options A  : $tree_inference_opts_a"
echo "Tree inf options B  : $tree_inference_opts_b"
echo "Date string         : $date_str"
echo "Tree desc A         : $tree_a_desc"
echo "Tree desc B         : $tree_b_desc"
echo "Tree A output name  : $tree_a_output_name"
echo "Tree B output name  : $tree_b_output_name"
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

if [ -f "$tree_a_output_name" ]; then
  echo "Tree $tree_a_output_name already exists"
  exit -1
fi

if [ -f "$tree_b_output_name" ]; then
  echo "Tree $tree_b_output_name already exists"
  exit -1
fi

if [ ! -f "$genomes_file_a" ]; then
  echo "$genomes_file_a does not exist"
  exit -1
fi

if [ ! -f "$genomes_file_b" ]; then
  echo "$genomes_file_b does not exist"
  exit -1
fi

if [ ! -f "$annotation_file" ]; then
  echo "$annotation_file does not exist"
  exit -1
fi

$CLEAN_FASTA -ml $minimal_length_a -r $char_ratio -t $seq_type $genomes_file_a $workdir/genomes_clean_a.fasta
$CLEAN_FASTA -ml $minimal_length_b -r $char_ratio -t $seq_type $genomes_file_b $workdir/genomes_clean_b.fasta

#only cluster if cutoot between 0 and 1:

$CDHIT -c $clustering_cutoff_a -d 0 -g 1 -i $workdir/genomes_clean_a.fasta -o $workdir/genomes_clean_a_cdhit.fasta
$CDHIT -c $clustering_cutoff_b -d 0 -g 1 -i $workdir/genomes_clean_b.fasta -o $workdir/genomes_clean_b_cdhit.fasta

mkdir $workdir/trees_a
mkdir $workdir/trees_b

cp $workdir/genomes_clean_a_cdhit.fasta $workdir/trees_a/
cp $workdir/genomes_clean_b_cdhit.fasta $workdir/trees_b/

$MK_TREE "$mafft_opts_a" "$msa_pro_opts_a" "-$tree_inference_opts_a" .fasta $workdir/trees_a/
$MK_TREE "$mafft_opts_b" "$msa_pro_opts_b" "-$tree_inference_opts_b" .fasta $workdir/trees_b/

$VIPR_X $workdir/trees_a/deco/genomes_clean_a_cdhit_ni_mafft_tree_raxml_d.xml $annotation_file $tree_a_output_name "$tree_a_desc"
$VIPR_X $workdir/trees_b/deco/genomes_clean_b_cdhit_ni_mafft_tree_raxml_d.xml $annotation_file $tree_b_output_name "$tree_b_desc"
